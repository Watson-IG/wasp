#!/usr/bin/env python3
"""
Cumulative-coverage plot per locus (no sci-notation anywhere).

Y-axis: % of locus bases ≥ depth  (ticks 0-100, axis spans 0-105)  
X-axis: depth (log scale) auto-ticks, plain numbers 1-1000+
USAGE (single line)
	python3 plot_depth_cdf.py --depth sample_per-base-depth.bed.gz --loci loci.bed --out depth_cdf.png
"""

import argparse, gzip
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ---------- helpers ----------
def load_loci(p: Path) -> pd.DataFrame:
	return (pd.read_csv(p, sep="\t", header=None,
	                    names=["chrom","start","end","name"],
	                    dtype={"chrom": str})
	          .sort_values(["chrom","start"])
	          .reset_index(drop=True))

def depths_per_locus(depth_gz: Path, loci_df: pd.DataFrame):
	by_chr = {}
	for idx, r in loci_df.iterrows():
		by_chr.setdefault(r.chrom, []).append((idx, r.start, r.end))
	depth_lists = [[] for _ in range(len(loci_df))]

	with gzip.open(depth_gz, "rt") as fh:
		for ln in fh:
			chrom, s, e, d = ln.strip().split("\t")
			s = int(s); e = int(e); d = int(d)
			for idx, L, R in by_chr.get(chrom, []):
				if e > L and s < R:                          # overlap
					depth_lists[idx].extend([d]*(min(e,R)-max(s,L)))

	for idx, lst in enumerate(depth_lists):
		yield loci_df.loc[idx, "name"], lst

def cdf(depths: list[int]):
	if not depths: return [], []
	cnt, total, rem = {}, len(depths), len(depths)
	for d in depths: cnt[d] = cnt.get(d, 0) + 1
	x, y = [], []
	for d in sorted(cnt):
		x.append(d);  y.append(100.0 * rem / total)
		rem -= cnt[d]
	return x, y

# ---------- main plot ----------
def plot(depth_gz: Path, loci_df: pd.DataFrame, out_png: Path):
	plt.figure(figsize=(5,4))
	for name, depths in depths_per_locus(depth_gz, loci_df):
		x, y = cdf(depths)
		if x: plt.plot(x, y, label=name, lw=1)

	ax = plt.gca()
	ax.set_xscale("log")

	# ---- disable sci-notation on X ----
	sf = ScalarFormatter()
	sf.set_scientific(False)
	sf.set_powerlimits((0, 0))
	ax.xaxis.set_major_formatter(sf)

	# ---- Y-axis: fixed 0-105, ticks every 25 ----
	ax.set_ylim(0, 105)
	ax.set_yticks([0, 25, 50, 75, 100])
	ax.set_yticklabels(["0", "25", "50", "75", "100"])

	ax.set_xlabel("HiFi Coverage (log)")
	ax.set_ylabel("% of the Locus")
	plt.legend(fontsize="xx-small", bbox_to_anchor=(1.02, 0.5),
	           loc="center left", borderaxespad=0)
	plt.tight_layout()
	plt.savefig(out_png, dpi=300)
	plt.close()

# ---------- CLI ----------
if __name__ == "__main__":
	p = argparse.ArgumentParser()
	p.add_argument("--depth", required=True)    # mosdepth per-base BED.gz
	p.add_argument("--loci",  required=True)    # 4-col BED
	p.add_argument("--out",   required=True)    # output PNG/PDF/SVG
	args = p.parse_args()
	plot(Path(args.depth), load_loci(Path(args.loci)), Path(args.out))
