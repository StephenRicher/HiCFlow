#!/usr/bin/env python3

"""Replace selected Nextflow work files with sparse files while preserving size and timestamps."""

import argparse
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "suffixes",
        nargs="+",
        help="Filename suffixes to sparsify, for example trimmed.fastq.gz truncated.fastq.gz.",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("analysis/work"),
        help="Nextflow work directory (default: analysis/work).",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("analysis/sparsified-work-files.tsv"),
        help="Metadata manifest path (default: analysis/sparsified-work-files.tsv).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List matched files without modifying them.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    work_dir = args.work_dir.resolve()

    if not work_dir.is_dir():
        raise SystemExit(f"Work directory does not exist: {work_dir}")

    suffixes = tuple(suffix.lstrip("*") for suffix in args.suffixes)
    files = sorted(
        path
        for path in work_dir.rglob("*")
        if path.is_file()
        and not path.is_symlink()
        and path.name.endswith(suffixes)
        and path.stat().st_size > 0
    )

    if not files:
        raise SystemExit(f"No non-empty files matched: {', '.join(suffixes)}")

    total_size = sum(path.stat().st_size for path in files)

    if args.dry_run:
        for path in files:
            print(f"{path}\t{path.stat().st_size}")
        print(f"Matched {len(files)} files ({total_size:,} logical bytes)")
        return

    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    with args.manifest.open("a", encoding="utf-8") as manifest:
        if args.manifest.stat().st_size == 0:
            manifest.write("path\tsize\tatime_ns\tmtime_ns\n")

        for path in files:
            stat = path.stat()
            manifest.write(
                f"{path}\t{stat.st_size}\t{stat.st_atime_ns}\t{stat.st_mtime_ns}\n"
            )

            with path.open("r+b") as handle:
                handle.truncate(0)
                handle.truncate(stat.st_size)

            os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
            print(f"Sparsified {path} ({stat.st_size:,} logical bytes)")

    print(f"Sparsified {len(files)} files ({total_size:,} logical bytes)")
    print(f"Metadata written to {args.manifest}")


if __name__ == "__main__":
    main()
