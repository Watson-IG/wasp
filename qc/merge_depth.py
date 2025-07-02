#!/usr/bin/env python3

import os
import gzip
import sys

def find_depth_files(results_dir='.'):
    """
    Crawl subdirectories of `results_dir` for:
    <sampleId>/stats/<sampleId>_ccs_to_ref-based_regions-depth.bed.gz
    Yields tuples of (sampleId, filepath).
    """
    for entry in os.scandir(results_dir):
        if entry.is_dir():
            sample_id = entry.name
            stats_dir = os.path.join(entry.path, 'stats')
            if os.path.isdir(stats_dir):
                for fname in os.listdir(stats_dir):
                    if fname.endswith('_ccs_to_ref-based_regions-depth.bed.gz'):
                        yield (sample_id, os.path.join(stats_dir, fname))

def parse_depth_file(filepath):
    """
    Parse a gzipped depth file and return a dict of {region: average_depth}.
    Expects file columns: chrom, start, end, region, depth.
    """
    region_depths = {}
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            # need at least 5 columns: chrom, start, end, region, depth
            if len(fields) < 5:
                continue

            region = fields[3]            # the region name
            try:
                depth = float(fields[4])  # the depth
            except ValueError:
                depth = 0.0

            region_depths[region] = depth

    return region_depths

def main():
    # Use command-line argument if provided; default to current working directory.
    results_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()

    all_depths = {}
    all_regions = set()

    # collect per-sample region→depth mappings
    for sample_id, depth_file in find_depth_files(results_dir):
        depth_dict = parse_depth_file(depth_file)
        all_depths[sample_id] = depth_dict
        all_regions.update(depth_dict.keys())

    # sort regions and build header
    sorted_regions = sorted(all_regions)
    header = ['sampleId'] + sorted_regions
    lines = ["\t".join(header)]

    # write one row per sample
    for sample_id in sorted(all_depths.keys()):
        row = [sample_id]
        for region in sorted_regions:
            # default to 0.0 if missing
            row.append(str(all_depths[sample_id].get(region, 0.0)))
        lines.append("\t".join(row))

    # write summary TSV
    output_path = os.path.join(results_dir, "depth_summary.tsv")
    with open(output_path, "w") as f:
        f.write("\n".join(lines))

if __name__ == '__main__':
    main()
