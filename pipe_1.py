import pandas as pd
import pysam
from pathlib import Path
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys


class Sample:

    def __init__(self, sample_id, r1, r2, run):
        self.id = str(sample_id)
        self.run = run

        self.r1 = Path(r1)
        self.r2 = Path(r2)

        self.base = Path("output") / run / self.id

        self.trimmed_dir = self.base / "trim"
        self.bam_dir = self.base / "bam"
        self.plots_dir = self.base / "plots"

        for d in [self.trimmed_dir, self.bam_dir, self.plots_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.r1_trim = self.trimmed_dir / "R1.fastq.gz"
        self.r2_trim = self.trimmed_dir / "R2.fastq.gz"

        self.bam = self.bam_dir / f"{self.id}.bam"

    def trim(self):
        subprocess.run([
            "fastp",
            "-i", str(self.r1),
            "-I", str(self.r2),
            "-o", str(self.r1_trim),
            "-O", str(self.r2_trim)
        ], check=True)

    def align(self, reference):

        cmd = (
            f"bwa mem {reference} {self.r1_trim} {self.r2_trim} "
            f"| samtools view -b - "
            f"| samtools sort -o {self.bam}"
        )

        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(["samtools", "index", str(self.bam)], check=True)

    def coverage(self, window=1_000_000):

        cov = {}

        with pysam.AlignmentFile(self.bam, "rb") as bam:
            for col in bam.pileup(min_mapping_quality=20):
                chrom = col.reference_name
                w = col.pos // window

                key = (chrom, w)
                cov[key] = cov.get(key, 0) + col.nsegments

        rows = []
        for (chrom, w), val in cov.items():
            rows.append([chrom, w * window, val / window])

        return pd.DataFrame(rows, columns=["chr", "start", "coverage"])


class Pipeline:

    def __init__(self, input_root, reference):

        self.input_root = Path(input_root)
        self.reference = reference

        self.samples = self.find_samples()

    def find_samples(self):

        samples = []

        for run_dir in self.input_root.iterdir():

            if not run_dir.is_dir():
                continue

            run_name = run_dir.name

            r1_files = list(run_dir.glob("*R1*.fastq*"))

            for r1 in r1_files:

                r2 = Path(str(r1).replace("R1", "R2"))

                if not r2.exists():
                    continue

                sample_id = r1.stem.split("_")[0]

                samples.append(Sample(sample_id, r1, r2, run_name))

        return samples

    def run_preprocessing(self):

        for s in self.samples:
            print(f"\nProcessing {s.run}/{s.id}")

            s.trim()
            s.align(self.reference)

    def compute_zscore(self, dfs):

        tables = []

        for sample_id, df in dfs.items():
            df = df.rename(columns={"coverage": str(sample_id)})
            tables.append(df)

        merged = tables[0]

        for t in tables[1:]:
            merged = merged.merge(t, on=["chr", "start"], how="inner")

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

    def add_genome_pos(self, df):

        df = df.copy()

        df["chr"] = df["chr"].astype(str).str.replace("chr", "")

        chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "M"]
        df["chr"] = pd.Categorical(df["chr"], categories=chrom_order, ordered=True)

        df = df.sort_values(["chr", "start"])

        chrom_sizes = df.groupby("chr")["start"].max().cumsum().shift(fill_value=0)

        df["genome_pos"] = df.apply(
            lambda row: row["start"] + chrom_sizes[row["chr"]],
            axis=1
        )

        return df

    def plot_zscores(self, z_df):

        out = Path("plots")
        out.mkdir(parents=True, exist_ok=True)

        df = self.add_genome_pos(z_df)

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

    def run(self):

        self.run_preprocessing()

        dfs = {}

        for s in self.samples:
            print(f"Coverage {s.id}")

            dfs[s.id] = s.coverage()

        z_df = self.compute_zscore(dfs)

        z_df.to_csv("zscore.tsv", sep="\t", index=False)

        print(z_df.head())

        self.plot_zscores(z_df)


if __name__ == "__main__":

    input_root = sys.argv[1]
    reference = sys.argv[2]

    pipeline = Pipeline(input_root, reference)
    pipeline.run()