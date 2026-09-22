#!/usr/bin/env python
"""Compare two jgi depth tables (v1.4.0 vs v1.5.0), ignoring column order.

Usage: compare_depths.py old.depth.txt new.depth.txt
"""
import sys
import pandas as pd

old = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
new = pd.read_csv(sys.argv[2], sep="\t", index_col=0)

if set(old.columns) != set(new.columns):
    print("Different columns:")
    print("  only in old:", sorted(set(old.columns) - set(new.columns)))
    print("  only in new:", sorted(set(new.columns) - set(old.columns)))
    sys.exit(1)
if set(old.index) != set(new.index):
    print(f"Different contigs: {len(set(old.index) ^ set(new.index))} not shared")
    sys.exit(1)

new = new.loc[old.index, old.columns]
diff = (old - new).abs().max()
print("Max absolute difference per column:")
print(diff.to_string())
print("\nOK" if (diff < 1e-3).all() else "\nDIFFERENCES FOUND")
