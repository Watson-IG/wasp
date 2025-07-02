#!/usr/bin/env python3
"""
Facet-plot mosdepth per-base depth across user-defined loci (4-column BED).

USAGE (one line, no escapes):
	python3 make_depth_plot.py --depth CDSR-01131-01_SH-5_ccs_to_ref-based_per-base-depth.bed.gz --loci new_IG_Loci.bed --out depth_facet.png
"""

import argparse
import gzip
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def load_loci(path: Path) -> pd.DataFrame:
	"""
	Read a 4-column BED (chrom, start, end, name) into a DataFrame
	with columns ['chrom','start','end','name'], sorted by chrom/start.
	"""
	df = pd.read_csv(
		path,
		sep="\t",
		header=None,
		names=["chrom", "start", "end", "name"],
		dtype={"chrom": str},
	)
	df.sort_values(["chrom", "start"], inplace=True)
	df.reset_index(drop=True, inplace=True)
	return df

def collect_depth(depth_gz: Path, loci_df: pd.DataFrame) -> list[list[tuple[int, int]]]:
	"""
	Scan the gzipped per-base depth BED and return:
	  facet_data[i] = list of (position, depth) tuples for locus i.
	"""
	# Initialize a list-of-lists for each locus index
	facet_data: list[list[tuple[int, int]]] = [[] for _ in range(len(loci_df))]
	# Build an index: chrom -> list of (locus_index, start, end)
	by_chrom: dict[str, list[tuple[int, int, int]]] = {}
	for idx, row in loci_df.iterrows():
		by_chrom.setdefault(row.chrom, []).append((idx, row.start, row.end))

	# Open the gzipped per-base file and stream through it line by line
	with gzip.open(depth_gz, "rt") as fh:
		for ln in fh:
			chrom, s_str, e_str, d_str = ln.strip().split("\t")
			s = int(s_str)
			e = int(e_str)
			d = int(d_str)
			# If this chromosome has any loci of interest
			for idx, L, R in by_chrom.get(chrom, []):
				# Check overlap: [s, e) overlaps [L, R) if e > L and s < R
				if e > L and s < R:
					facet_data[idx].append((s, d))
					facet_data[idx].append((e, d))

	# Sort and dedupe each locus's depth points
	for lst in facet_data:
		lst[:] = sorted(set(lst), key=lambda t: t[0])
	return facet_data

def plot_facets(
	depth_gz: Path,
	loci_df: pd.DataFrame,
	depth_data: list[list[tuple[int, int]]],
	out_png: Path,
) -> None:
	"""
	Given a DataFrame of loci and the per-locus depth_data, draw one panel per locus.
	Each panel’s x-axis will show multiple tick marks (no scientific notation).
	The y-axis is fixed from 0 to 300.
	"""
	n = len(loci_df)
	fig, axes = plt.subplots(
		nrows=n,
		ncols=1,
		figsize=(12, 1.5 * n),
		sharey=True,
		constrained_layout=True,
	)
	# If there's only one locus, axes is not a list by default—wrap it
	if n == 1:
		axes = [axes]

	for ax, (idx, row) in zip(axes, loci_df.iterrows()):
		pts = depth_data[idx]
		if pts:
			x_vals, y_vals = zip(*pts)
			ax.plot(x_vals, y_vals, linewidth=0.7)
		else:
			ax.text(0.5, 0.5, "no coverage", ha="center", va="center", transform=ax.transAxes)

		# Force the exact genomic range from start to end
		ax.set_xlim(row.start, row.end)

		# Turn off scientific notation on the x-axis
		ax.ticklabel_format(style="plain", axis="x")

		# Panel title = the 4th column (locus name)
		ax.set_title(row["name"], fontsize=9, pad=2)

		# X-axis label shows chrom:start-end
		ax.set_xlabel(f"{row.chrom}:{row.start}-{row.end}")

		# Y-axis label
		ax.set_ylabel("depth")

	# Fix y-axis range 0-300 for every panel
	for ax in axes:
        ax.set_ylim(0, 500)
        ax.set_yticks([0, 250, 500])

	# Overall title at the top
	plt.suptitle(depth_gz.name, fontsize=11)

	# Save to the specified file (PNG/PDF/SVG) and close
	plt.savefig(out_png, dpi=300)
	plt.close(fig)

def main():
	parser = argparse.ArgumentParser(description="Facet-plot mosdepth depth over loci.")
	parser.add_argument("--depth", required=True, help="mosdepth per-base BED.gz")
	parser.add_argument("--loci", required=True, help="4-column BED: chrom, start, end, name")
	parser.add_argument("--out", required=True, help="output image filename (e.g. PNG)")
	args = parser.parse_args()

	loci_df = load_loci(Path(args.loci))
	depth_data = collect_depth(Path(args.depth), loci_df)
	plot_facets(Path(args.depth), loci_df, depth_data, Path(args.out))

if __name__ == "__main__":
	main()
