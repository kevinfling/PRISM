#!/usr/bin/env python3
"""Summarize the fixed 20-replicate protocol; standard library only.

Intervals describe Monte Carlo variation under an independent-replicate model.
They are approximate t intervals for non-normal squared-error observations.
Nested E prefixes are summarized separately, never treated as independent groups.
"""
import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict


def summarize(path, output):
    groups = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            key = row["control"], int(row["E"])
            values = [float(row[name]) for name in
                      ("squared_error", "expected_squared_error", "variance_squared_error")]
            if any(not math.isfinite(x) or x < 0 for x in values):
                raise ValueError("nonfinite or negative squared-error moment")
            if (int(row["n"]), int(row["q"]), float(row["t"]), float(row["sigma"])) != (64, 3, 16, 0.25):
                raise ValueError("unexpected protocol parameters")
            if int(row["seed"]) != 20261021:
                raise ValueError("unexpected master seed")
            control_index = ("zero", "white", "fgn_h03", "fgn_h07").index(row["control"])
            replicate = int(row["replicate"])
            if int(row["stream"]) != 4 * replicate + control_index:
                raise ValueError("unexpected logical stream")
            groups[key].append((replicate, *values))
    expected_groups = {(c, e) for c in ("zero", "white", "fgn_h03", "fgn_h07")
                       for e in (1, 4, 16, 64, 256)}
    if set(groups) != expected_groups:
        raise ValueError("expected four controls at five ensemble counts")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("control", "E", "replicates", "mean_squared_error", "expected_squared_error",
                     "observed_expected_ratio", "empirical_standard_error", "approx_t95_lower",
                     "approx_t95_upper", "theoretical_standard_error"))
    for (control, count), rows in sorted(groups.items()):
        if len(rows) != 20 or {r[0] for r in rows} != set(range(20)):
            raise ValueError("expected exactly 20 distinct replicates per control/count")
        if len({(r[2], r[3]) for r in rows}) != 1:
            raise ValueError("theoretical moments differ across replicates")
        values = [r[1] for r in rows]
        mean = statistics.mean(values)
        se = statistics.stdev(values) / math.sqrt(20)
        expected, variance = rows[0][2:]
        writer.writerow((control, count, 20, mean, expected, mean / expected if expected else "",
                         se, mean - 2.0930240544 * se, mean + 2.0930240544 * se,
                         math.sqrt(variance / 20)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    args = parser.parse_args()
    try:
        summarize(args.csv_path, sys.stdout)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"invalid experiment input: {error}\n")


if __name__ == "__main__":
    main()
