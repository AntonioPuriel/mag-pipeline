#!/usr/bin/env bash
# Hand-checked test of bin/gene_abundance.py and bin/contig_profiles.py.
#
# Assemblies A and B have samples s1 and s2; assembly C only s1 (as with
# --assembly_mode per_sample). Catalog genes:
#   A|c1_1: members in A and B             -> mean across assemblies
#   A|c2_1: two members on A|c2            -> sum within the assembly
#   C|c9_1: members in C and A; s2 only in A -> s2 from A alone
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
bin="$here/../bin"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

printf "contigName\tcontigLen\ttotalAvgDepth\ts1.sorted.bam\ts1.sorted.bam-var\ts2.sorted.bam\ts2.sorted.bam-var\n" > A.depth.txt
printf "c1\t5000\t30\t10\t1\t20\t1\nc2\t3000\t6\t4\t1\t2\t1\n" >> A.depth.txt
printf "contigName\tcontigLen\ttotalAvgDepth\ts1.sorted.bam\ts1.sorted.bam-var\ts2.sorted.bam\ts2.sorted.bam-var\n" > B.depth.txt
printf "c1\t5000\t30\t12\t1\t18\t1\n" >> B.depth.txt
# C as written by other versions of jgi_summarize_bam_contig_depths: no .bam
printf "contigName\tcontigLen\ttotalAvgDepth\ts1.sorted\ts1.sorted-var\n" > C.depth.txt
printf "c9\t2000\t6\t6\t1\n" >> C.depth.txt

printf "A|c1_1\tA|c1_1\nA|c1_1\tB|c1_1\nA|c2_1\tA|c2_1\nA|c2_1\tA|c2_2\nC|c9_1\tC|c9_1\nC|c9_1\tA|c2_3\n" > clusters.tsv

python3 "$bin/gene_abundance.py" --clusters clusters.tsv --depths A.depth.txt B.depth.txt C.depth.txt \
    --out-depth depth.tsv.gz --out-tpm tpm.tsv.gz --out-stats stats.tsv

printf "## emapper\n#query\tseed_ortholog\tKEGG_ko\nA|c1_1\tx\tko:K00001,ko:K00002\nA|c2_1\tx\tko:K00001\nC|c9_1\tx\t-\n" > chunk001.emapper.annotations

python3 "$bin/contig_profiles.py" --tpm tpm.tsv.gz --eggnog chunk001.emapper.annotations --outdir .

python3 - <<'EOF'
import pandas as pd
d = pd.read_csv("depth.tsv.gz", sep="\t", index_col=0)
assert list(d.columns) == ["s1", "s2"], f"sample names: {list(d.columns)}"
expected = {"A|c1_1": (11, 19), "A|c2_1": (8, 4), "C|c9_1": (5, 2)}
for gene, (s1, s2) in expected.items():
    got = tuple(d.loc[gene, ["s1", "s2"]])
    assert got == (s1, s2), f"{gene}: depth {got}, expected {(s1, s2)}"

t = pd.read_csv("tpm.tsv.gz", sep="\t", index_col=0)
assert abs(t.loc["A|c1_1", "s1"] - 11 / 24 * 1e6) < 1, t
assert abs(t["s2"].sum() - 1e6) < 1, t

ko = pd.read_csv("ko_per_sample.tsv", sep="\t").set_index(["ko", "sample"])
assert ko.loc[("K00001", "s1"), "n_genes"] == 2
assert abs(ko.loc[("K00001", "s1"), "tpm"] - 19 / 24 * 1e6) < 1, ko
assert abs(ko.loc[("K00002", "s1"), "tpm"] - 11 / 24 * 1e6) < 1, ko
print("contig scripts: all checks passed")
EOF
