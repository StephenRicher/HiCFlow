#!/usr/bin/env python3

"""Convert GFF3 transcripts overlapping configured regions to BED12."""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple
from urllib.parse import unquote


TRANSCRIPT_FEATURES = {
    "transcript",
    "mRNA",
    "lnc_RNA",
    "ncRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "tRNA",
}

BIOTYPE_COLOURS = {
    "protein_coding": "0,0,180",
    "lncRNA": "180,0,0",
    "lincRNA": "180,0,0",
    "antisense": "180,0,0",
    "processed_transcript": "180,0,0",
    "pseudogene": "128,128,128",
}


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")

    return path.open(encoding="utf-8")


def parse_attributes(value: str) -> Dict[str, str]:
    attributes = {}

    for field in value.rstrip(";").split(";"):
        if not field:
            continue

        key, separator, raw_value = field.partition("=")
        if separator:
            attributes[key] = unquote(raw_value)

    return attributes


def read_regions(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    regions = defaultdict(list)

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected at least three columns."
                )

            chromosome, start, end = fields[:3]
            regions[chromosome].append((int(start), int(end)))

    if not regions:
        raise ValueError(f"No regions found in {path}.")

    return dict(regions)


def overlaps_regions(
    chromosome: str,
    start: int,
    end: int,
    regions: Dict[str, List[Tuple[int, int]]],
) -> bool:
    return any(start < region_end and end > region_start for region_start, region_end in regions.get(chromosome, []))


def attribute(
    attributes: Dict[str, str],
    *keys: str,
    default: Optional[str] = None,
) -> Optional[str]:
    for key in keys:
        value = attributes.get(key)
        if value:
            return value

    return default


def colour_for_biotype(biotype: str) -> str:
    if biotype in BIOTYPE_COLOURS:
        return BIOTYPE_COLOURS[biotype]

    if "pseudogene" in biotype:
        return BIOTYPE_COLOURS["pseudogene"]

    return "0,0,0"


def remove_transcript_prefix(value: str) -> str:
    prefix = "transcript:"
    return value[len(prefix):] if value.startswith(prefix) else value


def convert(gff3: Path, regions_path: Path, output: Path) -> None:
    regions = read_regions(regions_path)
    transcripts: Dict[str, Dict[str, object]] = {}
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open_text(gff3) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chromosome, _, feature, start_value, end_value, _, strand, _, attributes_value = fields

            if chromosome not in regions:
                continue

            start = int(start_value) - 1
            end = int(end_value)

            if not overlaps_regions(chromosome, start, end, regions):
                continue

            attributes = parse_attributes(attributes_value)

            if feature in TRANSCRIPT_FEATURES:
                transcript_id = attribute(attributes, "transcript_id", "ID")
                if not transcript_id:
                    continue

                transcript_id = remove_transcript_prefix(transcript_id)
                transcripts[transcript_id] = {
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "name": attribute(
                        attributes,
                        "gene_name",
                        "Name",
                        "gene_id",
                        default=transcript_id,
                    ),
                    "biotype": attribute(
                        attributes,
                        "gene_type",
                        "gene_biotype",
                        "transcript_type",
                        "transcript_biotype",
                        default="unknown",
                    ),
                }
                continue

            if feature != "exon":
                continue

            parents = attribute(attributes, "Parent")
            if not parents:
                transcript_id = attribute(attributes, "transcript_id")
                parents = transcript_id

            if not parents:
                continue

            for parent in parents.split(","):
                transcript_id = remove_transcript_prefix(parent)
                exons[transcript_id].append((start, end))

    with output.open("w", encoding="utf-8") as handle:
        for transcript_id in sorted(
            transcripts,
            key=lambda item: (
                str(transcripts[item]["chromosome"]),
                int(transcripts[item]["start"]),
                int(transcripts[item]["end"]),
                item,
            ),
        ):
            transcript = transcripts[transcript_id]
            transcript_exons = sorted(set(exons.get(transcript_id, [])))

            if not transcript_exons:
                transcript_exons = [
                    (
                        int(transcript["start"]),
                        int(transcript["end"]),
                    )
                ]

            chrom_start = min(start for start, _ in transcript_exons)
            chrom_end = max(end for _, end in transcript_exons)
            block_sizes = ",".join(str(end - start) for start, end in transcript_exons)
            block_starts = ",".join(str(start - chrom_start) for start, _ in transcript_exons)

            fields = [
                str(transcript["chromosome"]),
                str(chrom_start),
                str(chrom_end),
                str(transcript["name"]),
                "0",
                str(transcript["strand"]),
                str(chrom_start),
                str(chrom_end),
                colour_for_biotype(str(transcript["biotype"])),
                str(len(transcript_exons)),
                block_sizes,
                block_starts,
            ]
            handle.write("\t".join(fields) + "\n")

    if not output.stat().st_size:
        raise ValueError(
            f"No transcript annotations from {gff3} overlapped regions in {regions_path}."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gff3", type=Path)
    parser.add_argument("regions", type=Path)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    convert(args.gff3, args.regions, args.output)


if __name__ == "__main__":
    main()
