# WASP Toolkit

This repository provides scripts and utilities for immune receptor analysis. A Python
package exposes a `wasptk-read-support` command that computes read support metrics
from a mapped BAM and an annotation table.

## Installation

```bash
pip install -e .
```

## Usage

```bash
wasptk-read-support <allele_annotation.csv> <mapped.bam> <output.csv>
```

Optional flags allow overriding the column names in the annotation table using
`--contig-col`, `--start-col`, `--end-col`, and `--gene-col`.

### Output columns

For each row in the annotation table the command appends the following columns:

* `Total_Positions` – number of reference positions examined
* `Total_Reads_by_Positions` – sum of per-position coverage
* `Average_Coverage` – mean coverage across the region
* `Mismatched_Positions` – count of positions with >20% mismatches
* `Matched_Positions` – count of positions with >80% matches
* `Position_Mismatches` – colon separated mismatches at each position
* `Position_Matches` – colon separated matches at each position
* `Percent_Accuracy` – `(Matched_Positions / Total_Positions) * 100`
* `Positions_With_At_Least_10x_Coverage` – number of positions with coverage ≥10
* `Fully_Spanning_Reads` – reads that span the entire region
* `Fully_Spanning_Reads_100%_Match` – fully spanning reads with a perfect match

Constant-region genes (`IGKC`, `IGLC`, `TRAC`, `TRBC`, `TRDC`, `TRGC`, `IGHC`) include
additional coverage breakdowns by 10× thresholds and per‑exon match/spanning counts.
