import pysam
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import sys

def get_bams(input_dir):
    return list(Path(input_dir).rglob("*.bam"))

def compute_coverage(bam_path, window=1_000_000):
    cov = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for col in bam.pileup(min_mapping_quality=20):
            chrom = col.reference_name
            w = col.pos // window
            key = (chrom, w)
            cov[key] = cov.get(key, 0) + 1

    rows = []
    for (chrom, w), val in cov.items():
        rows.append([chrom, w * window, val])

    return pd.DataFrame(rows, columns=["chr", "start", "coverage"])

def zscore(dfs):
    tables = []

    for sample_id, df in dfs.items():
        df = df.rename(columns={"coverage": sample_id})
        tables.append(df)

    merged = tables[0]

    for t in tables[1:]:
        merged = merged.merge(t, on=["chr", "start"], how="outer")

    sample_cols = list(dfs.keys())

    coords = merged[["chr", "start"]].reset_index(drop=True)

    numeric = merged[sample_cols].fillna(0)

    values = numeric.to_numpy(dtype=float)

    mean = np.mean(values, axis=1, keepdims=True)
    std = np.std(values, axis=1, keepdims=True)

    std[std == 0] = 1

    z = (values - mean) / std

    z_df = pd.DataFrame(z, columns=[f"{c}_z" for c in sample_cols])

    return pd.concat([coords.reset_index(drop=True), z_df.reset_index(drop=True)], axis=1)

def add_genome_pos(df):

    df = df.copy()

    df["chr"] = df["chr"].str.replace("chr", "")

    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "M"]
    df["chr"] = pd.Categorical(df["chr"], categories=chrom_order, ordered=True)

    df = df.sort_values(["chr", "start"])

    chrom_sizes = df.groupby("chr")["start"].max().cumsum().shift(fill_value=0)

    df["genome_pos"] = df.apply(
        lambda row: row["start"] + chrom_sizes[row["chr"]],
        axis=1
    )
    return df

def plot(z_df):
    out = Path("plots")
    out.mkdir(exist_ok=True)

    df = add_genome_pos(z_df)

    for col in df.columns:
        if not col.endswith("_z"):
            continue

        plt.figure(figsize=(12, 4))
        plt.scatter(df["genome_pos"], df[col], s=2)

        plt.axhline(0)
        plt.axhline(2, linestyle="--")
        plt.axhline(-2, linestyle="--")
        plt.xlabel("Genome position")
        plt.ylabel("Z-score")

        out_file = out / f"{col}.png"

        if out_file.exists():
            out_file.unlink()

        plt.savefig(out_file, dpi=150)
        plt.close()

def main():
    input_dir = sys.argv[1]

    bams = get_bams(input_dir)

    dfs = {}

    for bam in bams:
        print(bam)
        dfs[bam.stem] = compute_coverage(bam)

    z_df = zscore(dfs)

    z_df.to_csv("zscore.tsv", sep="\t", index=False)

    plot(z_df)

if __name__ == "__main__":
    main()