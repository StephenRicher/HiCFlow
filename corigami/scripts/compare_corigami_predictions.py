#!/usr/bin/env python3

"""Compare C.Origami wild-type and deletion predictions."""

import argparse
import json
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


def load_prediction(path: Path) -> np.ndarray:
    """Load and validate a C.Origami prediction matrix."""

    matrix = np.squeeze(np.load(path))

    if matrix.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional prediction in {path}, "
            f"found shape {matrix.shape}."
        )

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"Expected a square prediction matrix in {path}, "
            f"found shape {matrix.shape}."
        )

    if not np.isfinite(matrix).all():
        raise ValueError(f"Prediction contains non-finite values: {path}")

    return matrix.astype(np.float64, copy=False)


def percentile_limit(
    values: np.ndarray,
    percentile: float,
    minimum: float = 1e-12,
) -> float:
    """Return a positive percentile-based plotting limit."""

    finite = np.abs(values[np.isfinite(values)])
    finite = finite[finite > 0]

    if finite.size == 0:
        return minimum

    return max(float(np.percentile(finite, percentile)), minimum)


def format_coordinate(value: float) -> str:
    """Format genomic coordinates for plot axes."""

    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.2f} Mb"

    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f} kb"

    return f"{value:.0f}"


def plot_matrix(
    matrix: np.ndarray,
    output: Path,
    title: str,
    chromosome: str,
    window_start: int,
    window_end: int,
    colourmap: str,
    vmin: float,
    vmax: float,
    deletion: Optional[Tuple[int, int]],
    deletion_label: Optional[str],
    colourbar_label: str,
    dpi: int,
) -> None:
    """Plot a prediction or difference matrix."""

    figure, axis = plt.subplots(figsize=(8, 7))

    image = axis.imshow(
        matrix,
        cmap=colourmap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
        interpolation="nearest",
        extent=[
            window_start,
            window_end,
            window_end,
            window_start,
        ],
        aspect="equal",
    )

    ticks = np.linspace(window_start, window_end, 5)
    labels = [format_coordinate(value) for value in ticks]

    axis.set_xticks(ticks)
    axis.set_xticklabels(labels)
    axis.set_yticks(ticks)
    axis.set_yticklabels(labels)
    axis.set_xlabel(f"{chromosome} coordinate (hg38)")
    axis.set_ylabel(f"{chromosome} coordinate (hg38)")
    axis.set_title(title, loc="left")

    if deletion is not None:
        if not deletion_label:
            raise ValueError(
                "A deletion label is required when displaying a deletion."
            )

        deletion_start, deletion_end = deletion
        deletion_midpoint = (deletion_start + deletion_end) / 2

        axis.annotate(
            deletion_label,
            xy=(deletion_midpoint, 1),
            xycoords=axis.get_xaxis_transform(),
            xytext=(0, 24),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="black",
            arrowprops={
                "arrowstyle": "-|>",
                "color": "black",
                "linewidth": 1.8,
                "mutation_scale": 12,
            },
            annotation_clip=False,
        )

    colourbar = figure.colorbar(
        image,
        ax=axis,
        fraction=0.046,
        pad=0.04,
    )
    colourbar.set_label(colourbar_label)

    figure.tight_layout(rect=[0, 0, 1, 0.94])
    figure.savefig(
        output,
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(figure)


def extract_viewpoint_profile(
    matrix: np.ndarray,
    window_start: int,
    window_size: int,
    viewpoint_start: int,
    viewpoint_end: int,
) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    """Extract a matrix profile from the bin containing a viewpoint."""

    matrix_size = matrix.shape[0]
    bin_size = window_size // matrix_size

    if bin_size * matrix_size != window_size:
        raise ValueError(
            f"Window size {window_size} is not divisible by "
            f"matrix size {matrix_size}."
        )

    viewpoint_midpoint = (viewpoint_start + viewpoint_end) // 2
    viewpoint_index = (viewpoint_midpoint - window_start) // bin_size

    if viewpoint_index < 0 or viewpoint_index >= matrix_size:
        raise ValueError(
            f"Viewpoint {viewpoint_start}-{viewpoint_end} falls outside "
            f"the prediction window beginning at {window_start}."
        )

    bin_starts = (
        window_start
        + np.arange(matrix_size, dtype=np.int64) * bin_size
    )
    bin_midpoints = bin_starts + (bin_size / 2)

    # C.Origami matrices should be symmetric. Averaging the corresponding
    # row and column makes the viewpoint robust to small numerical asymmetry.
    profile = (
        matrix[viewpoint_index, :]
        + matrix[:, viewpoint_index]
    ) / 2

    viewpoint_bin_start = int(bin_starts[viewpoint_index])
    viewpoint_bin_end = viewpoint_bin_start + bin_size

    distances = bin_midpoints - (
        viewpoint_bin_start + viewpoint_bin_end
    ) / 2

    return (
        distances,
        profile,
        viewpoint_index,
        viewpoint_bin_start,
        viewpoint_bin_end,
    )


def write_viewpoint_table(
    output: Path,
    chromosome: str,
    window_start: int,
    window_size: int,
    distances: np.ndarray,
    wild_type_profile: np.ndarray,
    deletion_profile: np.ndarray,
) -> None:
    """Write WT, deletion and difference viewpoint profiles."""

    matrix_size = wild_type_profile.size
    bin_size = window_size // matrix_size

    bin_starts = (
        window_start
        + np.arange(matrix_size, dtype=np.int64) * bin_size
    )
    bin_ends = bin_starts + bin_size
    difference = deletion_profile - wild_type_profile

    with output.open("w", encoding="utf-8") as handle:
        handle.write(
            "chromosome\tstart\tend\tdistance\t"
            "wild_type\tdeletion\tdeletion_minus_WT\n"
        )

        for values in zip(
            bin_starts,
            bin_ends,
            distances,
            wild_type_profile,
            deletion_profile,
            difference,
        ):
            start, end, distance, wt, edited, change = values
            handle.write(
                f"{chromosome}\t{start}\t{end}\t"
                f"{distance:.1f}\t{wt:.10g}\t"
                f"{edited:.10g}\t{change:.10g}\n"
            )


def plot_viewpoint_comparison(
    wild_type: np.ndarray,
    deletion: np.ndarray,
    output: Path,
    title: str,
    chromosome: str,
    window_start: int,
    window_size: int,
    deletion_start: int,
    deletion_end: int,
    dpi: int,
) -> dict[str, int]:
    """Plot WT, deletion and directional difference viewpoint profiles."""

    (
        distances,
        wild_type_profile,
        viewpoint_index,
        viewpoint_bin_start,
        viewpoint_bin_end,
    ) = extract_viewpoint_profile(
        matrix=wild_type,
        window_start=window_start,
        window_size=window_size,
        viewpoint_start=deletion_start,
        viewpoint_end=deletion_end,
    )

    (
        deletion_distances,
        deletion_profile,
        deletion_viewpoint_index,
        deletion_bin_start,
        deletion_bin_end,
    ) = extract_viewpoint_profile(
        matrix=deletion,
        window_start=window_start,
        window_size=window_size,
        viewpoint_start=deletion_start,
        viewpoint_end=deletion_end,
    )

    if not np.array_equal(distances, deletion_distances):
        raise ValueError("WT and deletion viewpoint distances differ.")

    if viewpoint_index != deletion_viewpoint_index:
        raise ValueError("WT and deletion viewpoint bins differ.")

    difference = deletion_profile - wild_type_profile

    figure, (profile_axis, difference_axis) = plt.subplots(
        2,
        1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )

    profile_axis.plot(
        distances,
        wild_type_profile,
        label="Predicted WT",
        color="#4C72B0",
        linewidth=1.5,
    )
    profile_axis.plot(
        distances,
        deletion_profile,
        label="Predicted deletion",
        color="#C44E52",
        linewidth=1.5,
    )
    profile_axis.axvline(
        0,
        color="black",
        linewidth=0.8,
        alpha=0.5,
    )
    profile_axis.set_ylabel("C.Origami predicted contact score")
    profile_axis.legend(loc="upper right")

    difference_axis.plot(
        distances,
        difference,
        color="black",
        linewidth=1.2,
    )
    difference_axis.fill_between(
        distances,
        0,
        difference,
        where=difference >= 0,
        color="firebrick",
        alpha=0.4,
        interpolate=True,
        label="Stronger after deletion",
    )
    difference_axis.fill_between(
        distances,
        0,
        difference,
        where=difference < 0,
        color="steelblue",
        alpha=0.4,
        interpolate=True,
        label="Weaker after deletion",
    )
    difference_axis.axhline(
        0,
        color="black",
        linewidth=0.8,
    )
    difference_axis.axvline(
        0,
        color="black",
        linewidth=0.8,
        alpha=0.5,
    )
    difference_axis.set_ylabel("Deletion − WT")
    difference_axis.legend(loc="upper right")

    difference_axis.set_xlabel(
        "Distance from B1 deletion viewpoint (bp, hg38)\n"
        f"{chromosome}:{viewpoint_bin_start}-{viewpoint_bin_end}"
    )

    figure.suptitle(
        title,
        x=0.01,
        ha="left",
    )
    figure.tight_layout()
    figure.savefig(
        output,
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(figure)

    return {
        "viewpoint_index": viewpoint_index,
        "viewpoint_bin_start": viewpoint_bin_start,
        "viewpoint_bin_end": viewpoint_bin_end,
    }


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wild-type", type=Path, required=True)
    parser.add_argument("--deletion", type=Path, required=True)
    parser.add_argument("--cell-line", required=True)
    parser.add_argument("--clone", required=True)
    parser.add_argument("--allele", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--window-start", type=int, required=True)
    parser.add_argument("--window-size", type=int, default=2_097_152)
    parser.add_argument("--deletion-start", type=int, required=True)
    parser.add_argument("--deletion-end", type=int, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--deletion-label", required=True, help="Label displayed above the deletion-minus-WT matrix.")
    parser.add_argument("--prediction-percentile", type=float, default=99.5)
    parser.add_argument("--difference-percentile", type=float, default=99.0)
    parser.add_argument(
        "--difference-limit",
        type=float,
        help=(
            "Fixed symmetric difference colour limit. If omitted, the limit "
            "is calculated from the requested percentile."
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)

    return parser.parse_args()


def main() -> None:
    """Create prediction and difference matrices and plots."""

    args = parse_args()

    wild_type = load_prediction(args.wild_type)
    deletion = load_prediction(args.deletion)

    if wild_type.shape != deletion.shape:
        raise ValueError(
            "Prediction shapes do not match: "
            f"WT {wild_type.shape}, deletion {deletion.shape}."
        )

    difference = deletion - wild_type
    window_end = args.window_start + args.window_size

    combined_predictions = np.concatenate(
        [
            wild_type.ravel(),
            deletion.ravel(),
        ]
    )
    finite_predictions = combined_predictions[
        np.isfinite(combined_predictions)
    ]

    prediction_tail = (100 - args.prediction_percentile) / 2
    prediction_vmin = float(
        np.percentile(
            finite_predictions,
            prediction_tail,
        )
    )
    prediction_vmax = float(
        np.percentile(
            finite_predictions,
            100 - prediction_tail,
        )
    )

    if prediction_vmin >= prediction_vmax:
        raise ValueError(
            "Prediction plotting limits are invalid: "
            f"{prediction_vmin} >= {prediction_vmax}."
        )

    difference_limit = (
        args.difference_limit
        if args.difference_limit is not None
        else percentile_limit(
            difference,
            args.difference_percentile,
        )
    )

    if difference_limit <= 0:
        raise ValueError(
            f"Difference limit must be positive, found {difference_limit}."
        )

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    np.save(f"{prefix}-wild-type.npy", wild_type)
    np.save(f"{prefix}-deletion.npy", deletion)
    np.save(f"{prefix}-deletion-minus-WT.npy", difference)

    deletion_interval = (
        args.deletion_start,
        args.deletion_end,
    )
    identifier = (
        f"{args.cell_line} clone {args.clone}, "
        f"allele {args.allele}"
    )
    coordinates = (
        f"{args.chromosome}:{args.window_start}-{window_end}"
    )

    plot_matrix(
        matrix=wild_type,
        output=Path(f"{prefix}-wild-type.svg"),
        title=f"{args.cell_line} predicted wild type\n{coordinates}",
        chromosome=args.chromosome,
        window_start=args.window_start,
        window_end=window_end,
        colourmap="Purples",
        vmin=prediction_vmin,
        vmax=prediction_vmax,
        deletion=None,
        deletion_label=None,
        colourbar_label="C.Origami predicted contact score",
        dpi=args.dpi,
    )

    plot_matrix(
        matrix=deletion,
        output=Path(f"{prefix}-deletion.svg"),
        title=f"{identifier}, predicted deletion\n{coordinates}",
        chromosome=args.chromosome,
        window_start=args.window_start,
        window_end=window_end,
        colourmap="Purples",
        vmin=prediction_vmin,
        vmax=prediction_vmax,
        deletion=None,
        deletion_label=None,
        colourbar_label="C.Origami predicted contact score",
        dpi=args.dpi,
    )

    plot_matrix(
        matrix=difference,
        output=Path(f"{prefix}-deletion-minus-WT.svg"),
        title=f"{identifier}, predicted deletion minus WT\n{coordinates}",
        chromosome=args.chromosome,
        window_start=args.window_start,
        window_end=window_end,
        colourmap="RdBu_r",
        vmin=-difference_limit,
        vmax=difference_limit,
        deletion=deletion_interval,
        deletion_label=args.deletion_label,
        colourbar_label=(
            "Deletion − WT predicted contact change\n"
            "Blue: weaker after deletion; red: stronger after deletion"
        ),
        dpi=args.dpi,
    )

    (
        viewpoint_distances,
        wild_type_viewpoint,
        viewpoint_index,
        viewpoint_bin_start,
        viewpoint_bin_end,
    ) = extract_viewpoint_profile(
        matrix=wild_type,
        window_start=args.window_start,
        window_size=args.window_size,
        viewpoint_start=args.deletion_start,
        viewpoint_end=args.deletion_end,
    )

    (
        deletion_viewpoint_distances,
        deletion_viewpoint,
        deletion_viewpoint_index,
        deletion_viewpoint_bin_start,
        deletion_viewpoint_bin_end,
    ) = extract_viewpoint_profile(
        matrix=deletion,
        window_start=args.window_start,
        window_size=args.window_size,
        viewpoint_start=args.deletion_start,
        viewpoint_end=args.deletion_end,
    )

    if not np.array_equal(
        viewpoint_distances,
        deletion_viewpoint_distances,
    ):
        raise ValueError("WT and deletion viewpoint coordinates differ.")

    if viewpoint_index != deletion_viewpoint_index:
        raise ValueError("WT and deletion viewpoints use different bins.")

    write_viewpoint_table(
        output=Path(f"{prefix}-viewpoint.tsv"),
        chromosome=args.chromosome,
        window_start=args.window_start,
        window_size=args.window_size,
        distances=viewpoint_distances,
        wild_type_profile=wild_type_viewpoint,
        deletion_profile=deletion_viewpoint,
    )

    plot_viewpoint_comparison(
        wild_type=wild_type,
        deletion=deletion,
        output=Path(f"{prefix}-viewpoint.svg"),
        title=(
            f"{identifier}, predicted B1 viewpoint\n"
            f"{coordinates}"
        ),
        chromosome=args.chromosome,
        window_start=args.window_start,
        window_size=args.window_size,
        deletion_start=args.deletion_start,
        deletion_end=args.deletion_end,
        dpi=args.dpi,
    )

    metadata = {
        "cell_line": args.cell_line,
        "clone": args.clone,
        "allele": args.allele,
        "chromosome": args.chromosome,
        "assembly": "hg38",
        "window_start": args.window_start,
        "window_end": window_end,
        "window_size": args.window_size,
        "matrix_shape": list(wild_type.shape),
        "bin_size": args.window_size / wild_type.shape[0],
        "deletion_start": args.deletion_start,
        "deletion_end": args.deletion_end,
        "deletion_width": args.deletion_end - args.deletion_start,
        "operation": "deletion - wild_type",
        "positive_values": "stronger predicted contact after deletion",
        "negative_values": "weaker predicted contact after deletion",
        "prediction_vmin": prediction_vmin,
        "prediction_vmax": prediction_vmax,
        "difference_colour_limit": difference_limit,
        "prediction_percentile": args.prediction_percentile,
        "difference_percentile": args.difference_percentile,
        "viewpoint_type": "C.Origami prediction bin overlapping deletion midpoint",
        "viewpoint_index": viewpoint_index,
        "viewpoint_bin_start": viewpoint_bin_start,
        "viewpoint_bin_end": viewpoint_bin_end,
        "viewpoint_bin_size": args.window_size // wild_type.shape[0],
        "viewpoint_profile": "mean of matrix row and corresponding column",
    }

    Path(f"{prefix}-metadata.json").write_text(
        json.dumps(
            metadata,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
