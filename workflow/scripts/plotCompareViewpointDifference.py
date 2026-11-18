#!/usr/bin/env python3

"""Plot two adjusted viewpoint profiles and their directional difference."""

import argparse
import sys
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd


COLUMNS = [
    "view_chromosome",
    "view_start",
    "view_end",
    "chromosome",
    "start",
    "end",
    "score",
]

COORDINATE_COLUMNS = COLUMNS[:-1]


def parse_viewpoint(value: str) -> Tuple[str, int, int]:
    """Parse chromosome:start-end viewpoint coordinates."""

    chromosome, coordinates = value.split(":", 1)
    start, end = coordinates.split("-", 1)

    return chromosome, int(start) - 1, int(end)


def read_profile(path: Path, score_column: str) -> pd.DataFrame:
    """Read and aggregate a hicPlotViewpoint bedGraph."""

    data = pd.read_csv(
        path,
        sep="\t",
        names=COLUMNS,
        dtype={
            "view_chromosome": str,
            "view_start": int,
            "view_end": int,
            "chromosome": str,
            "start": int,
            "end": int,
            "score": float,
        },
    )

    if data.empty:
        raise ValueError(f"Viewpoint profile is empty: {path}")

    return (
        data.groupby(
            ["chromosome", "start", "end"],
            as_index=False,
            sort=True,
        )
        .agg(
            {
                "view_chromosome": "first",
                "view_start": "min",
                "view_end": "max",
                "score": "sum",
            }
        )
        .rename(columns={"score": score_column})
    )

def combine_profiles(
    group1_path: Path,
    group2_path: Path,
    group1: str,
    group2: str,
    viewpoint: Tuple[str, int, int],
) -> pd.DataFrame:
    """Join profiles by genomic bin and calculate group2 minus group1."""

    group1_data = read_profile(group1_path, group1)
    group2_data = read_profile(group2_path, group2)

    combined = group1_data.merge(
        group2_data,
        on=COORDINATE_COLUMNS,
        how="outer",
        validate="one_to_one",
    )

    missing_group1 = combined[group1].isna().sum()
    missing_group2 = combined[group2].isna().sum()

    if missing_group1 or missing_group2:
        raise ValueError(
            "The two profiles do not contain identical genomic bins: "
            f"{missing_group1} missing from {group1}, "
            f"{missing_group2} missing from {group2}."
        )

    combined["midpoint"] = (combined["start"] + combined["end"]) / 2
    _, viewpoint_start, viewpoint_end = viewpoint
    combined["viewpoint_midpoint"] = (viewpoint_start + viewpoint_end) / 2
    combined["distance"] = (
        combined["midpoint"] - combined["viewpoint_midpoint"]
    )
    combined["difference"] = combined[group2] - combined[group1]

    return combined


def plot_comparison(
    data: pd.DataFrame,
    group1: str,
    group2: str,
    viewpoint: Tuple[str, int, int],
    build: Optional[str],
    title: Optional[str],
    out: Path,
    dpi: int,
) -> None:
    """Plot adjusted profiles and their directional difference."""

    figure, (profile_axis, difference_axis) = plt.subplots(
        2,
        1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )

    profile_axis.plot(
        data["distance"],
        data[group1],
        label=group1,
        linewidth=1.5,
    )
    profile_axis.plot(
        data["distance"],
        data[group2],
        label=group2,
        linewidth=1.5,
    )
    profile_axis.axvline(0, color="black", linestyle="--", alpha=0.5)
    profile_axis.set_ylabel("Adjusted interactions")
    profile_axis.legend(loc="upper right")

    difference_axis.plot(
        data["distance"],
        data["difference"],
        color="black",
        linewidth=1.5,
    )
    difference_axis.fill_between(
        data["distance"],
        0,
        data["difference"],
        where=data["difference"] >= 0,
        color="firebrick",
        alpha=0.35,
        interpolate=True,
        label=f"Higher in {group2}",
    )
    difference_axis.fill_between(
        data["distance"],
        0,
        data["difference"],
        where=data["difference"] < 0,
        color="steelblue",
        alpha=0.35,
        interpolate=True,
        label=f"Higher in {group1}",
    )
    difference_axis.axhline(0, color="black", linewidth=0.8)
    difference_axis.axvline(
        0,
        color="black",
        linestyle="--",
        alpha=0.5,
    )
    difference_axis.set_ylabel(f"{group2} - {group1}")
    difference_axis.legend(loc="upper right")

    chromosome, start, end = viewpoint
    viewpoint_label = f"{chromosome}:{start + 1}-{end}"

    xlabel = f"Distance from viewpoint at {viewpoint_label} (bp)"
    if build:
        xlabel += f" ({build})"

    difference_axis.set_xlabel(xlabel)

    if title:
        figure.suptitle(title, x=0.01, ha="left")

    figure.tight_layout()
    figure.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("group1_bedgraph", type=Path)
    parser.add_argument("group2_bedgraph", type=Path)
    parser.add_argument("--viewpoint", required=True, help="Viewpoint coordinates in chromosome:start-end format.")
    parser.add_argument("--group1", required=True)
    parser.add_argument("--group2", required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--build")
    parser.add_argument("--title")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    viewpoint = parse_viewpoint(args.viewpoint)

    data = combine_profiles(
        args.group1_bedgraph,
        args.group2_bedgraph,
        args.group1,
        args.group2,
        viewpoint,
    )

    output_columns = [
        "chromosome",
        "start",
        "end",
        args.group1,
        args.group2,
        "difference",
    ]

    data[output_columns].to_csv(
        args.table,
        sep="\t",
        index=False,
    )

    plot_comparison(
        data=data,
        group1=args.group1,
        group2=args.group2,
        viewpoint=viewpoint,
        build=args.build,
        title=args.title,
        out=args.out,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    sys.exit(main())
