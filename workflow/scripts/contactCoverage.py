#!/usr/bin/env python3

"""Plot cumulative contact coverage for Hi-C matrices."""

import argparse
import sys
from collections.abc import Sequence

import matplotlib.pyplot as plt
from hicmatrix import HiCMatrix as hm

from utilities import createMainParent, setDefaults

__version__ = "2.0.0"


def plot_coverage(
    matrices: Sequence[str],
    resolutions: Sequence[int],
    sample: str,
    out: str,
    n_bins: int,
    dpi: int,
    font_size: float,
) -> None:
    """Plot the proportion of bins exceeding each contact count."""

    if len(matrices) != len(resolutions):
        raise ValueError(
            f"Received {len(matrices)} matrices but {len(resolutions)} resolutions."
        )

    if not matrices:
        raise ValueError("At least one matrix is required.")

    plt.rcParams.update({"font.size": font_size})
    figure, axis = plt.subplots(figsize=(16, 8))

    x_limit = 0.0

    for matrix, resolution in zip(matrices, resolutions, strict=True):
        contacts_per_bin = hm.hiCMatrix(matrix).matrix.sum(axis=1).A1

        cumulative_counts, bin_edges, _ = axis.hist(
            contacts_per_bin,
            n_bins,
            density=True,
            cumulative=-1,
            histtype="step",
            label=str(resolution),
        )

        below_threshold = bin_edges[:-1][cumulative_counts < 0.1]
        cutoff = below_threshold.min() if below_threshold.size else bin_edges.max()
        x_limit = max(x_limit, float(cutoff))

    axis.set_xlim(0, x_limit)
    axis.set_title(sample, loc="left")
    axis.legend(title="Bin size", loc="upper right")
    axis.set_xlabel("Contacts per bin")
    axis.set_ylabel("Proportion")

    figure.tight_layout()
    figure.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def parse_args():
    epilog = "Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)"
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=epilog,
        parents=[createMainParent(verbose=False, version=__version__)],
    )

    parser.add_argument(
        "matrices",
        nargs="+",
        help="Hi-C matrices in H5 format.",
    )
    parser.add_argument(
        "--resolutions",
        nargs="+",
        type=int,
        required=True,
        help="Matrix resolutions in the same order as the matrix arguments.",
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample or group name displayed in the plot.",
    )
    parser.add_argument(
        "--nBins",
        dest="n_bins",
        type=int,
        default=10_000,
        help="Number of histogram bins.",
    )
    parser.add_argument(
        "--fontSize",
        dest="font_size",
        type=float,
        default=12,
        help="Plot font size.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output resolution.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output plot path.",
    )
    parser.set_defaults(function=plot_coverage)

    return setDefaults(parser)


if __name__ == "__main__":
    args, function = parse_args()
    sys.exit(function(**vars(args)))
