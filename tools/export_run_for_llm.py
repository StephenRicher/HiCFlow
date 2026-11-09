#!/usr/bin/env python3

"""Export a HiCFlow run as a structured, LLM-friendly analysis package."""

import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


TEXT_SUFFIXES = {
    ".awk",
    ".bed",
    ".bed12",
    ".bedgraph",
    ".config",
    ".csv",
    ".gff",
    ".gff3",
    ".groovy",
    ".ini",
    ".json",
    ".links",
    ".log",
    ".md",
    ".nf",
    ".out",
    ".py",
    ".r",
    ".sh",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}

PLOT_SUFFIXES = {
    ".jpg",
    ".jpeg",
    ".pdf",
    ".png",
    ".svg",
    ".webp",
}

LARGE_DATA_SUFFIXES = {
    ".bam",
    ".bai",
    ".bt2",
    ".bw",
    ".fastq",
    ".fq",
    ".h5",
    ".hic",
    ".pkl",
    ".sam",
}

COMPRESSED_SUFFIXES = {
    ".gz",
    ".bgz",
    ".bz2",
    ".xz",
    ".zip",
}

SOURCE_PATTERNS = (
    "main.nf",
    "Dockerfile",
    "environment.yml",
    "environment-linux-64.lock",
    ".dockerignore",
    ".gitignore",
    "example/example.config",
    "modules/**/*.nf",
    "workflow/scripts/**/*.py",
    "workflow/scripts/**/*.R",
    "workflow/scripts/**/*.r",
    "workflow/scripts/**/*.sh",
    "workflow/scripts/**/*.awk",
    "tests/**/*.sh",
    "tests/**/*.py",
)

CONFIG_PATTERNS = (
    "*.config",
    "config/**/*",
)

OPTIONAL_SECTIONS = {
    "cscore": (
        "results/cscore",
        "results/Cscore",
        "cscore",
        "Cscore",
    ),
    "loops": (
        "results/loops",
        "loops",
    ),
}

SECRET_PATTERNS = (
    re.compile(
        r"(?i)\b(password|passwd|secret|token|api[_-]?key|access[_-]?key)"
        r"\s*[:=]\s*(['\"]?)([^'\"\s]+)\2"
    ),
    re.compile(r"(?i)\b(ghp|github_pat|glpat)-[A-Za-z0-9_-]+"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
)

LOG_ERROR_PATTERN = re.compile(
    r"(?i)\b("
    r"error|exception|traceback|failed|failure|killed|out of memory|oom|"
    r"segmentation fault|cannot|missing output|no such file"
    r")\b"
)

LOG_WARNING_PATTERN = re.compile(r"(?i)\b(warn|warning|deprecated)\b")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run",
        type=Path,
        help=(
            "Run directory, for example 'analysis'. It should contain "
            "results/, and may contain work/, config/, and analysis.config."
        ),
    )
    parser.add_argument(
        "--results",
        type=Path,
        help="Results directory. Defaults to RUN/results.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "Export directory. Defaults to "
            "RUN/llm_export_<UTC timestamp>."
        ),
    )
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path.cwd(),
        help="HiCFlow repository root. Defaults to the current directory.",
    )
    parser.add_argument(
        "--max-text-bytes",
        type=int,
        default=2_000_000,
        help="Maximum bytes copied from one text file. Default: 2000000.",
    )
    parser.add_argument(
        "--log-tail-lines",
        type=int,
        default=500,
        help="Number of trailing lines retained per log. Default: 500.",
    )
    parser.add_argument(
        "--copy-html",
        action="store_true",
        help="Copy HTML reports as well as plots and text files.",
    )
    parser.add_argument(
        "--hash-limit-bytes",
        type=int,
        default=100_000_000,
        help=(
            "Calculate SHA-256 only for files at or below this size. "
            "Default: 100000000."
        ),
    )
    return parser.parse_args()


def utc_timestamp() -> str:
    return datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def relative_to_or_absolute(path: Path, parent: Path) -> str:
    try:
        return str(path.relative_to(parent))
    except ValueError:
        return str(path)


def sorted_files(root: Path) -> List[Path]:
    if not root.exists():
        return []

    return sorted(
        (
            path
            for path in root.rglob("*")
            if path.is_file() and not path.is_symlink()
        ),
        key=lambda path: str(path),
    )


def format_bytes(value: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    amount = float(value)

    for unit in units:
        if amount < 1024 or unit == units[-1]:
            return "{:.2f} {}".format(amount, unit)
        amount /= 1024

    return "{} B".format(value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def redact_secrets(text: str) -> str:
    redacted = text

    for pattern in SECRET_PATTERNS:
        if pattern.groups >= 3:
            redacted = pattern.sub(
                lambda match: "{}=<REDACTED>".format(match.group(1)),
                redacted,
            )
        else:
            redacted = pattern.sub("<REDACTED>", redacted)

    return redacted


def read_text(path: Path, max_bytes: int) -> Tuple[str, bool]:
    size = path.stat().st_size
    truncated = size > max_bytes

    with path.open("rb") as handle:
        data = handle.read(max_bytes)

    return data.decode("utf-8", errors="replace"), truncated


def read_gzip_text(path: Path, max_bytes: int) -> Tuple[str, bool]:
    chunks: List[bytes] = []
    total = 0
    truncated = False

    with gzip.open(path, "rb") as handle:
        while total < max_bytes:
            chunk = handle.read(min(1024 * 1024, max_bytes - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)

        if handle.read(1):
            truncated = True

    return b"".join(chunks).decode("utf-8", errors="replace"), truncated


def tail_lines(path: Path, limit: int) -> List[str]:
    if limit <= 0:
        return []

    lines: List[str] = []

    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                lines.append(line.rstrip("\n"))
                if len(lines) > limit:
                    lines.pop(0)
    except OSError as exc:
        return ["Could not read log: {}".format(exc)]

    return lines


def run_command(
    command: Sequence[str],
    cwd: Optional[Path] = None,
) -> Dict[str, Any]:
    try:
        completed = subprocess.run(
            list(command),
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError as exc:
        return {
            "command": list(command),
            "returncode": None,
            "stdout": "",
            "stderr": str(exc),
        }

    return {
        "command": list(command),
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def git_metadata(repo: Path) -> Dict[str, Any]:
    commands = {
        "commit": ["git", "rev-parse", "HEAD"],
        "branch": ["git", "branch", "--show-current"],
        "status": ["git", "status", "--short"],
        "remote": ["git", "remote", "-v"],
        "describe": ["git", "describe", "--always", "--dirty", "--tags"],
        "diff_stat": ["git", "diff", "--stat"],
        "staged_diff_stat": ["git", "diff", "--cached", "--stat"],
    }

    metadata: Dict[str, Any] = {}

    for name, command in commands.items():
        result = run_command(command, cwd=repo)
        metadata[name] = {
            "returncode": result["returncode"],
            "value": result["stdout"].strip(),
            "error": result["stderr"].strip(),
        }

    return metadata


def write_json(path: Path, value: Any) -> None:
    ensure_directory(path.parent)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_text(path: Path, value: str) -> None:
    ensure_directory(path.parent)
    path.write_text(value, encoding="utf-8")


def copy_file(source: Path, destination: Path) -> None:
    ensure_directory(destination.parent)
    shutil.copy2(str(source), str(destination))


def find_source_files(repo: Path) -> List[Path]:
    files = set()

    for pattern in SOURCE_PATTERNS:
        for path in repo.glob(pattern):
            if path.is_file():
                files.add(path.resolve())

    return sorted(files, key=lambda path: str(path))


def export_source(repo: Path, output: Path, max_bytes: int) -> List[str]:
    source_root = output / "workflow_source"
    files = find_source_files(repo)
    manifest: List[str] = []
    combined_parts: List[str] = []

    for source in files:
        relative = source.relative_to(repo.resolve())
        destination = source_root / relative
        text, truncated = read_text(source, max_bytes)
        text = redact_secrets(text)

        if truncated:
            text += (
                "\n\n# EXPORT NOTE: file truncated at {} bytes.\n".format(
                    max_bytes
                )
            )

        write_text(destination, text)
        manifest.append(str(relative))

        combined_parts.extend(
            [
                "=" * 100,
                "FILE: {}".format(relative),
                "=" * 100,
                text.rstrip(),
                "",
            ]
        )

    write_text(
        output / "workflow_source_combined.txt",
        "\n".join(combined_parts).rstrip() + "\n",
    )
    write_text(
        output / "workflow_source_manifest.txt",
        "\n".join(manifest) + "\n",
    )

    return manifest


def export_configs(
    run_root: Path,
    output: Path,
    max_bytes: int,
) -> List[str]:
    exported = set()
    destination_root = output / "configuration"

    candidates = [
        run_root / "analysis.config",
        run_root / "nextflow.config",
    ]

    for pattern in CONFIG_PATTERNS:
        candidates.extend(run_root.glob(pattern))

    for source in sorted(set(candidates), key=lambda path: str(path)):
        if not source.is_file():
            continue

        relative = source.relative_to(run_root)
        destination = destination_root / relative
        text, truncated = read_text(source, max_bytes)
        text = redact_secrets(text)

        if truncated:
            text += "\n\n# EXPORT NOTE: file truncated.\n"

        write_text(destination, text)
        exported.add(str(relative))

    write_text(
        output / "configuration_manifest.txt",
        "\n".join(sorted(exported)) + "\n",
    )

    return sorted(exported)


def inventory_entry(
    path: Path,
    base: Path,
    hash_limit: int,
) -> Dict[str, Any]:
    stat = path.stat()
    suffixes = "".join(path.suffixes).lower()

    entry = {
        "path": relative_to_or_absolute(path, base),
        "size_bytes": stat.st_size,
        "size_human": format_bytes(stat.st_size),
        "modified_utc": datetime.datetime.utcfromtimestamp(
            stat.st_mtime
        ).isoformat() + "Z",
        "suffix": suffixes,
    }

    if stat.st_size <= hash_limit:
        entry["sha256"] = sha256(path)
    else:
        entry["sha256"] = None
        entry["hash_note"] = "Skipped because file exceeds hash limit."

    return entry


def write_inventory(
    results: Path,
    output: Path,
    hash_limit: int,
) -> List[Dict[str, Any]]:
    files = sorted_files(results)
    entries = [
        inventory_entry(path, results, hash_limit)
        for path in files
    ]

    write_json(output / "file_inventory.json", entries)

    with (output / "file_inventory.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "path",
                "size_bytes",
                "size_human",
                "modified_utc",
                "suffix",
                "sha256",
            ],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(entries)

    return entries


def section_from_path(relative: str) -> str:
    parts = Path(relative).parts

    if not parts:
        return "root"

    if len(parts) == 1:
        return parts[0]

    return "{}/{}".format(parts[0], parts[1])


def summarise_inventory(
    entries: Sequence[Dict[str, Any]],
    output: Path,
) -> None:
    extension_counts = Counter()
    section_counts = Counter()
    section_sizes = Counter()
    total_size = 0

    for entry in entries:
        suffix = entry["suffix"] or "<none>"
        section = section_from_path(entry["path"])
        size = int(entry["size_bytes"])

        extension_counts[suffix] += 1
        section_counts[section] += 1
        section_sizes[section] += size
        total_size += size

    lines = [
        "# Results inventory summary",
        "",
        "Total files: {}".format(len(entries)),
        "Total size: {}".format(format_bytes(total_size)),
        "",
        "## By results section",
        "",
        "Section\tFiles\tSize",
    ]

    for section in sorted(section_counts):
        lines.append(
            "{}\t{}\t{}".format(
                section,
                section_counts[section],
                format_bytes(section_sizes[section]),
            )
        )

    lines.extend(
        [
            "",
            "## By file type",
            "",
            "Suffix\tFiles",
        ]
    )

    for suffix, count in sorted(extension_counts.items()):
        lines.append("{}\t{}".format(suffix, count))

    write_text(output / "inventory_summary.txt", "\n".join(lines) + "\n")


def copy_plots(results: Path, output: Path) -> List[str]:
    destination_root = output / "plots"
    copied: List[str] = []

    for source in sorted_files(results):
        if source.suffix.lower() not in PLOT_SUFFIXES:
            continue

        relative = source.relative_to(results)
        copy_file(source, destination_root / relative)
        copied.append(str(relative))

    write_text(
        output / "plots_manifest.txt",
        "\n".join(copied) + "\n",
    )

    return copied


def copy_html_reports(results: Path, output: Path) -> List[str]:
    copied: List[str] = []
    destination_root = output / "html_reports"

    for source in sorted_files(results):
        if source.suffix.lower() != ".html":
            continue

        relative = source.relative_to(results)
        copy_file(source, destination_root / relative)
        copied.append(str(relative))

    write_text(
        output / "html_reports_manifest.txt",
        "\n".join(copied) + "\n",
    )

    return copied


def is_useful_text_file(path: Path) -> bool:
    suffix = path.suffix.lower()

    if suffix in TEXT_SUFFIXES:
        return True

    if suffix == ".gz":
        stem_suffix = Path(path.stem).suffix.lower()
        return stem_suffix in TEXT_SUFFIXES

    return False


def copy_tabular_outputs(
    results: Path,
    output: Path,
    max_bytes: int,
) -> List[Dict[str, Any]]:
    copied: List[Dict[str, Any]] = []
    destination_root = output / "tabular_outputs"

    for source in sorted_files(results):
        relative = source.relative_to(results)

        if "logs" in relative.parts:
            continue

        if not is_useful_text_file(source):
            continue

        try:
            if source.suffix.lower() == ".gz":
                text, truncated = read_gzip_text(source, max_bytes)
            else:
                text, truncated = read_text(source, max_bytes)
        except (OSError, EOFError) as exc:
            copied.append(
                {
                    "path": str(relative),
                    "status": "unreadable",
                    "error": str(exc),
                }
            )
            continue

        destination_name = relative
        if source.suffix.lower() == ".gz":
            destination_name = relative.with_suffix("")

        destination = destination_root / destination_name
        write_text(destination, redact_secrets(text))

        copied.append(
            {
                "path": str(relative),
                "exported_as": str(destination.relative_to(output)),
                "truncated": truncated,
                "original_size_bytes": source.stat().st_size,
            }
        )

    write_json(output / "tabular_outputs_manifest.json", copied)
    return copied


def parse_delimited_counts(path: Path) -> Optional[Dict[str, Any]]:
    if path.suffix.lower() not in {".bed", ".bedgraph", ".csv", ".tsv", ".links"}:
        return None

    line_count = 0
    data_line_count = 0
    chromosome_counts = Counter()

    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line_count += 1
                stripped = line.strip()

                if not stripped or stripped.startswith(("#", "track", "browser")):
                    continue

                data_line_count += 1
                fields = re.split(r"\t|,", stripped)
                if fields:
                    chromosome_counts[fields[0]] += 1
    except OSError:
        return None

    return {
        "path": str(path),
        "line_count": line_count,
        "data_line_count": data_line_count,
        "chromosome_counts": dict(sorted(chromosome_counts.items())),
    }


def summarise_structural_outputs(
    results: Path,
    output: Path,
) -> Dict[str, Any]:
    categories = {
        "tads": results / "tads",
        "loops": results / "loops",
        "cscore": results / "cscore",
        "comparisons": results / "matrix" / "compare",
        "viewpoints": results / "viewpoints",
    }
    summary: Dict[str, Any] = {}

    for name, root in categories.items():
        files = sorted_files(root)
        records = []

        for path in files:
            relative = path.relative_to(results)
            record = {
                "path": str(relative),
                "size_bytes": path.stat().st_size,
                "empty": path.stat().st_size == 0,
            }

            counts = parse_delimited_counts(path)
            if counts:
                record.update(
                    {
                        "line_count": counts["line_count"],
                        "data_line_count": counts["data_line_count"],
                        "chromosome_counts": counts["chromosome_counts"],
                    }
                )

            records.append(record)

        summary[name] = {
            "present": root.exists(),
            "file_count": len(files),
            "non_empty_file_count": sum(
                1 for path in files if path.stat().st_size > 0
            ),
            "files": records,
        }

    write_json(output / "structural_outputs_summary.json", summary)

    lines = ["# Structural outputs summary", ""]

    for name in ("tads", "loops", "cscore", "comparisons", "viewpoints"):
        value = summary[name]
        lines.extend(
            [
                "## {}".format(name),
                "",
                "Present: {}".format(value["present"]),
                "Files: {}".format(value["file_count"]),
                "Non-empty files: {}".format(value["non_empty_file_count"]),
                "",
            ]
        )

        if not value["files"]:
            lines.append("No outputs found.")
            lines.append("")
            continue

        lines.append("Path\tBytes\tData lines\tEmpty")
        for record in value["files"]:
            lines.append(
                "{}\t{}\t{}\t{}".format(
                    record["path"],
                    record["size_bytes"],
                    record.get("data_line_count", ""),
                    record["empty"],
                )
            )
        lines.append("")

    write_text(
        output / "structural_outputs_summary.txt",
        "\n".join(lines) + "\n",
    )

    return summary


def parse_key_value_file(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}

    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                fields = stripped.split("\t", 1)
                if len(fields) != 2:
                    fields = stripped.split(None, 1)

                if len(fields) == 2:
                    values[fields[0]] = fields[1]
    except OSError:
        pass

    return values


def summarise_comparison_metadata(
    results: Path,
    output: Path,
) -> List[Dict[str, Any]]:
    records = []

    for path in sorted(results.glob("matrix/compare/**/*-subtraction.tsv")):
        if not path.is_file():
            continue

        records.append(
            {
                "path": str(path.relative_to(results)),
                "metadata": parse_key_value_file(path),
            }
        )

    write_json(output / "comparison_metadata.json", records)

    lines = ["# Matrix comparison metadata", ""]

    if not records:
        lines.append("No subtraction metadata found.")
    else:
        for record in records:
            lines.append("## {}".format(record["path"]))
            lines.append("")
            for key, value in record["metadata"].items():
                lines.append("{}: {}".format(key, value))
            lines.append("")

    write_text(
        output / "comparison_metadata.txt",
        "\n".join(lines) + "\n",
    )

    return records


def read_trace(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error):
        return []


def summarise_pipeline_trace(
    results: Path,
    output: Path,
) -> Dict[str, Any]:
    candidates = sorted(
        results.glob("pipeline_info/*trace*.txt"),
        key=lambda path: path.stat().st_mtime,
    )

    if not candidates:
        summary = {
            "present": False,
            "message": "No Nextflow trace file found.",
        }
        write_json(output / "pipeline_trace_summary.json", summary)
        write_text(
            output / "pipeline_trace_summary.txt",
            "No Nextflow trace file found.\n",
        )
        return summary

    trace_path = candidates[-1]
    records = read_trace(trace_path)
    statuses = Counter()
    processes = Counter()
    cached = 0

    for record in records:
        status = record.get("status", "<unknown>")
        process = record.get("process", "<unknown>")
        statuses[status] += 1
        processes[process] += 1

        if record.get("cached", "").lower() in {"true", "yes", "1"}:
            cached += 1

    summary = {
        "present": True,
        "path": str(trace_path.relative_to(results)),
        "task_count": len(records),
        "cached_task_count": cached,
        "statuses": dict(sorted(statuses.items())),
        "processes": dict(sorted(processes.items())),
        "records": records,
    }

    write_json(output / "pipeline_trace_summary.json", summary)

    lines = [
        "# Nextflow trace summary",
        "",
        "Trace: {}".format(summary["path"]),
        "Tasks: {}".format(summary["task_count"]),
        "Cached tasks: {}".format(summary["cached_task_count"]),
        "",
        "## Statuses",
        "",
    ]

    for status, count in sorted(statuses.items()):
        lines.append("{}\t{}".format(status, count))

    lines.extend(["", "## Processes", ""])

    for process, count in sorted(processes.items()):
        lines.append("{}\t{}".format(process, count))

    write_text(
        output / "pipeline_trace_summary.txt",
        "\n".join(lines) + "\n",
    )

    return summary


def export_logs(
    results: Path,
    output: Path,
    tail_limit: int,
) -> Dict[str, Any]:
    log_root = results / "logs"
    destination_root = output / "logs"
    records = []
    error_lines = []
    warning_lines = []

    for source in sorted_files(log_root):
        if source.suffix.lower() != ".log":
            continue

        relative = source.relative_to(results)
        lines = tail_lines(source, tail_limit)
        text = "\n".join(lines) + ("\n" if lines else "")
        text = redact_secrets(text)
        destination = destination_root / relative
        write_text(destination, text)

        errors = [
            line for line in lines
            if LOG_ERROR_PATTERN.search(line)
        ]
        warnings = [
            line for line in lines
            if LOG_WARNING_PATTERN.search(line)
        ]

        records.append(
            {
                "path": str(relative),
                "original_size_bytes": source.stat().st_size,
                "exported_lines": len(lines),
                "possible_error_lines": len(errors),
                "possible_warning_lines": len(warnings),
            }
        )

        for line in errors:
            error_lines.append("{}\t{}".format(relative, line))

        for line in warnings:
            warning_lines.append("{}\t{}".format(relative, line))

    write_json(output / "logs_summary.json", records)
    write_text(
        output / "possible_log_errors.txt",
        "\n".join(error_lines) + ("\n" if error_lines else ""),
    )
    write_text(
        output / "possible_log_warnings.txt",
        "\n".join(warning_lines) + ("\n" if warning_lines else ""),
    )

    return {
        "file_count": len(records),
        "possible_error_line_count": len(error_lines),
        "possible_warning_line_count": len(warning_lines),
    }


def collect_qc_files(results: Path, output: Path) -> Dict[str, Any]:
    qc_root = results / "qc"
    files = sorted_files(qc_root)
    extension_counts = Counter()
    section_counts = Counter()
    total_size = 0

    records = []

    for path in files:
        relative = path.relative_to(results)
        suffix = "".join(path.suffixes).lower() or "<none>"
        section = relative.parts[1] if len(relative.parts) > 1 else "root"
        size = path.stat().st_size

        extension_counts[suffix] += 1
        section_counts[section] += 1
        total_size += size

        records.append(
            {
                "path": str(relative),
                "size_bytes": size,
                "empty": size == 0,
            }
        )

    summary = {
        "present": qc_root.exists(),
        "file_count": len(files),
        "total_size_bytes": total_size,
        "by_section": dict(sorted(section_counts.items())),
        "by_suffix": dict(sorted(extension_counts.items())),
        "files": records,
    }

    write_json(output / "qc_summary.json", summary)

    lines = [
        "# QC summary",
        "",
        "Present: {}".format(summary["present"]),
        "Files: {}".format(summary["file_count"]),
        "Size: {}".format(format_bytes(summary["total_size_bytes"])),
        "",
        "## Sections",
        "",
    ]

    for section, count in sorted(section_counts.items()):
        lines.append("{}\t{}".format(section, count))

    lines.extend(["", "## Files", "", "Path\tBytes\tEmpty"])

    for record in records:
        lines.append(
            "{}\t{}\t{}".format(
                record["path"],
                record["size_bytes"],
                record["empty"],
            )
        )

    write_text(output / "qc_summary.txt", "\n".join(lines) + "\n")
    return summary


def detect_optional_sections(
    run_root: Path,
    results: Path,
    output: Path,
) -> Dict[str, Any]:
    status = {}

    for name, candidates in OPTIONAL_SECTIONS.items():
        found = []

        for candidate in candidates:
            result_candidate = results / candidate
            run_candidate = run_root / candidate

            for path in (result_candidate, run_candidate):
                if path.exists():
                    found.append(
                        {
                            "path": relative_to_or_absolute(path, run_root),
                            "file_count": len(sorted_files(path)),
                            "non_empty_file_count": sum(
                                1
                                for file_path in sorted_files(path)
                                if file_path.stat().st_size > 0
                            ),
                        }
                    )

        status[name] = {
            "present": bool(found),
            "locations": found,
            "note": (
                "Optional output found."
                if found
                else "Optional output was not produced or was disabled."
            ),
        }

    write_json(output / "optional_outputs.json", status)

    lines = ["# Optional outputs", ""]

    for name, value in sorted(status.items()):
        lines.extend(
            [
                "## {}".format(name),
                "",
                "Present: {}".format(value["present"]),
                value["note"],
                "",
            ]
        )

        for location in value["locations"]:
            lines.append(
                "{}: {} files, {} non-empty".format(
                    location["path"],
                    location["file_count"],
                    location["non_empty_file_count"],
                )
            )
        lines.append("")

    write_text(output / "optional_outputs.txt", "\n".join(lines) + "\n")
    return status


def matrix_manifest(results: Path, output: Path) -> List[Dict[str, Any]]:
    records = []

    for path in sorted(results.glob("matrix/**/*.h5")):
        if not path.is_file():
            continue

        relative = path.relative_to(results)
        name = path.name

        matrix_type = "unknown"
        for candidate in (
            "LOESSdiff-medianFilter",
            "LOESSdiff-noFilter",
            "obs-exp",
            "obsExp",
            "adjIF1",
            "adjIF2",
            "normalised",
            "KR",
            "raw",
            "zscore",
        ):
            if candidate.lower() in name.lower():
                matrix_type = candidate
                break

        records.append(
            {
                "path": str(relative),
                "filename": name,
                "matrix_type": matrix_type,
                "size_bytes": path.stat().st_size,
                "size_human": format_bytes(path.stat().st_size),
            }
        )

    write_json(output / "matrix_manifest.json", records)

    lines = [
        "# Matrix manifest",
        "",
        "Path\tType\tSize",
    ]

    for record in records:
        lines.append(
            "{}\t{}\t{}".format(
                record["path"],
                record["matrix_type"],
                record["size_human"],
            )
        )

    write_text(output / "matrix_manifest.txt", "\n".join(lines) + "\n")
    return records


def write_package_readme(
    output: Path,
    run_root: Path,
    results: Path,
    source_count: int,
    config_count: int,
    inventory: Sequence[Dict[str, Any]],
    plots: Sequence[str],
    optional_status: Dict[str, Any],
) -> None:
    lines = [
        "# HiCFlow LLM analysis export",
        "",
        "Generated: {}".format(
            datetime.datetime.utcnow().isoformat() + "Z"
        ),
        "Run directory: `{}`".format(run_root),
        "Results directory: `{}`".format(results),
        "",
        "This package contains compact text summaries, selected result files, "
        "plots, logs, workflow source, configuration, and provenance for "
        "LLM-assisted analysis.",
        "",
        "It deliberately does not copy large FASTQ, BAM, H5, BigWig, index, "
        "or work-directory files. Those files are listed in the inventory.",
        "",
        "## Start here",
        "",
        "1. `run_overview.txt`",
        "2. `configuration_manifest.txt` and `configuration/`",
        "3. `workflow_source_combined.txt`",
        "4. `pipeline_trace_summary.txt`",
        "5. `comparison_metadata.txt`",
        "6. `structural_outputs_summary.txt`",
        "7. `qc_summary.txt`",
        "8. `possible_log_errors.txt`",
        "9. `plots_manifest.txt` and `plots/`",
        "",
        "## Package counts",
        "",
        "- Result files inventoried: {}".format(len(inventory)),
        "- Workflow/config source files: {}".format(source_count),
        "- Run configuration files: {}".format(config_count),
        "- Plots copied: {}".format(len(plots)),
        "",
        "## Optional analyses",
        "",
    ]

    for name, value in sorted(optional_status.items()):
        lines.append(
            "- {}: {}".format(
                name,
                "present" if value["present"] else "not present",
            )
        )

    lines.extend(
        [
            "",
            "The absence of Cscore or loop outputs may mean the feature was "
            "disabled, unsuitable for the analysed region, failed, or produced "
            "no reportable features. Consult the configuration and logs before "
            "drawing a biological conclusion.",
            "",
        ]
    )

    write_text(output / "README.md", "\n".join(lines))


def write_run_overview(
    output: Path,
    run_root: Path,
    results: Path,
    git: Dict[str, Any],
    matrix_records: Sequence[Dict[str, Any]],
    qc_summary: Dict[str, Any],
    structural_summary: Dict[str, Any],
    log_summary: Dict[str, Any],
    plot_count: int,
) -> None:
    commit = git.get("commit", {}).get("value", "")
    branch = git.get("branch", {}).get("value", "")
    describe = git.get("describe", {}).get("value", "")
    status = git.get("status", {}).get("value", "")

    lines = [
        "# HiCFlow run overview",
        "",
        "Run directory: {}".format(run_root),
        "Results directory: {}".format(results),
        "Generated: {}".format(
            datetime.datetime.utcnow().isoformat() + "Z"
        ),
        "",
        "## Workflow provenance",
        "",
        "Git commit: {}".format(commit or "<unavailable>"),
        "Git branch: {}".format(branch or "<unavailable>"),
        "Git description: {}".format(describe or "<unavailable>"),
        "Repository dirty: {}".format(bool(status)),
        "",
        "## Result counts",
        "",
        "Matrices: {}".format(len(matrix_records)),
        "Plots: {}".format(plot_count),
        "QC files: {}".format(qc_summary.get("file_count", 0)),
        "TAD files: {}".format(
            structural_summary.get("tads", {}).get("file_count", 0)
        ),
        "Loop files: {}".format(
            structural_summary.get("loops", {}).get("file_count", 0)
        ),
        "Cscore files: {}".format(
            structural_summary.get("cscore", {}).get("file_count", 0)
        ),
        "Comparison files: {}".format(
            structural_summary.get("comparisons", {}).get("file_count", 0)
        ),
        "Viewpoint files: {}".format(
            structural_summary.get("viewpoints", {}).get("file_count", 0)
        ),
        "",
        "## Log scan",
        "",
        "Logs exported: {}".format(log_summary.get("file_count", 0)),
        "Possible error lines: {}".format(
            log_summary.get("possible_error_line_count", 0)
        ),
        "Possible warning lines: {}".format(
            log_summary.get("possible_warning_line_count", 0)
        ),
        "",
        "Log matches are heuristic. Inspect the corresponding exported log "
        "before treating a matched line as a real failure.",
        "",
    ]

    write_text(output / "run_overview.txt", "\n".join(lines))


def main() -> int:
    args = parse_args()

    repo = args.repo.expanduser().resolve()
    run_root = args.run.expanduser().resolve()
    results = (
        args.results.expanduser().resolve()
        if args.results
        else run_root / "results"
    )
    output = (
        args.output.expanduser().resolve()
        if args.output
        else run_root / "llm_export_{}".format(utc_timestamp())
    )

    if not repo.is_dir():
        raise SystemExit("Repository directory does not exist: {}".format(repo))

    if not run_root.is_dir():
        raise SystemExit("Run directory does not exist: {}".format(run_root))

    if not results.is_dir():
        raise SystemExit("Results directory does not exist: {}".format(results))

    if output.exists() and any(output.iterdir()):
        raise SystemExit(
            "Output directory already exists and is not empty: {}".format(
                output
            )
        )

    ensure_directory(output)

    git = git_metadata(repo)
    write_json(output / "git_metadata.json", git)

    source_manifest = export_source(
        repo,
        output,
        args.max_text_bytes,
    )
    config_manifest = export_configs(
        run_root,
        output,
        args.max_text_bytes,
    )

    inventory = write_inventory(
        results,
        output,
        args.hash_limit_bytes,
    )
    summarise_inventory(inventory, output)

    plots = copy_plots(results, output)

    if args.copy_html:
        copy_html_reports(results, output)

    copy_tabular_outputs(
        results,
        output,
        args.max_text_bytes,
    )

    matrix_records = matrix_manifest(results, output)
    structural_summary = summarise_structural_outputs(results, output)
    comparison_metadata = summarise_comparison_metadata(results, output)
    trace_summary = summarise_pipeline_trace(results, output)
    qc_summary = collect_qc_files(results, output)
    log_summary = export_logs(
        results,
        output,
        args.log_tail_lines,
    )
    optional_status = detect_optional_sections(
        run_root,
        results,
        output,
    )

    write_run_overview(
        output,
        run_root,
        results,
        git,
        matrix_records,
        qc_summary,
        structural_summary,
        log_summary,
        len(plots),
    )

    write_package_readme(
        output,
        run_root,
        results,
        len(source_manifest),
        len(config_manifest),
        inventory,
        plots,
        optional_status,
    )

    export_metadata = {
        "generated_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "repo": str(repo),
        "run": str(run_root),
        "results": str(results),
        "output": str(output),
        "result_file_count": len(inventory),
        "plot_count": len(plots),
        "source_file_count": len(source_manifest),
        "configuration_file_count": len(config_manifest),
        "comparison_metadata_count": len(comparison_metadata),
        "trace_present": trace_summary.get("present", False),
        "optional_outputs": optional_status,
    }
    write_json(output / "export_metadata.json", export_metadata)

    print("Created HiCFlow LLM export: {}".format(output))
    print("Start with: {}".format(output / "README.md"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
