import pandas as pd
import pysam
from pathlib import Path
import numpy as np
import subprocess



output_dir = Path("output")
class Sample:
    
    def __init__(self, sample_id, r1_path, r2_path):
        self.id = sample_id
        self.r1 = Path(r1_path)
        self.r2 = Path(r2_path)
        self.sample_dir = output_dir / self.id
        
        self.trimmed_dir = Path("trimmed")
        self.bam_dir = Path("bam")
        self.plots_dir = Path("plots")

        for d in [self.trimmed_dir, self.bam_dir, self.plots_dir]:
            d.mkdir(parents=True, exist_ok=True)
            
        self.r1_trim = self.trimmed_dir / f"{self.id}_R1.trim.fastq.gz"
        self.r2_trim = self.trimmed_dir / f"{self.id}_R2.trim.fastq.gz"
        self.bam_file = self.bam_dir / f"{self.id}.sorted.bam"

    def trim(self):

        cmd = [
            "fastp",
            "-i", str(self.r1),
            "-I", str(self.r2),
            "-o", str(self.r1_trim),
            "-O", str(self.r2_trim),
            "--qualified_quality_phred", "20",
            "--length_required", "30"
        ]

        subprocess.run(cmd, check=True)

    def align(self, reference):

        cmd = (
            f"bwa mem {reference} {self.r1_trim} {self.r2_trim} "
            f"| samtools view -b - "
            f"| samtools sort -o {self.bam_file}"
        )

        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(["samtools", "index", str(self.bam_file)], check=True)

    def coverage_table(self, window_size=1_000_000):

        coverage = {}

        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for col in bam.pileup():
                chrom = col.reference_name
                window = col.pos // window_size
                coverage[(chrom, window)] = coverage.get((chrom, window), 0) + col.nsegments

        rows = []
        for (chrom, window), depth in coverage.items():
            start = window * window_size
            rows.append([chrom, start, depth])

        df = pd.DataFrame(rows, columns=["chr", "start", "coverage"])

        return df


class Pipeline:

    def __init__(self, samples_table, reference):
        self.samples_table = samples_table
        self.reference = reference
        self.samples = self.load_samples()

    def load_samples(self):
        df = pd.read_csv(self.samples_table, sep="\t")
        return [
            Sample(row["ID"], row["R1"], row["R2"])
            for _, row in df.iterrows()
        ]

    def run_preprocessing(self):
        for sample in self.samples:
            print(f"\n=== Processing {sample.id} ===")

            sample.trim()
            sample.align(self.reference)

    def compute_zscore(self, dfs):

        tables = []

        for sample_id, df in dfs.items():
            df = df[["chr", "start", "coverage"]].copy()
            df = df.rename(columns={"coverage": sample_id})
            tables.append(df)

        merged = tables[0]
        for t in tables[1:]:
            merged = merged.merge(t, on=["chr", "start"], how="inner")

        sample_cols = list(dfs.keys())

        coords = merged[["chr", "start"]].reset_index(drop=True)

        numeric = merged[sample_cols].apply(pd.to_numeric, errors="coerce")

        mask = numeric.notnull().all(axis=1)
        numeric = numeric[mask].reset_index(drop=True)
        coords = coords[mask].reset_index(drop=True)

        nonzero_mask = (numeric != 0).any(axis=1)
        numeric = numeric[nonzero_mask].reset_index(drop=True)
        coords = coords[nonzero_mask].reset_index(drop=True)

        values = numeric.to_numpy(dtype=float)
        mean = np.mean(values, axis=1, keepdims=True)
        std = np.std(values, axis=1, keepdims=True)

        std[std == 0] = 1

        z = (values - mean) / std

        z_cols = [f"{col}_z" for col in sample_cols]
        z_df = pd.DataFrame(z, columns=z_cols)

        result = pd.concat([coords, z_df], axis=1)

        return result

    def run(self):

        self.run_preprocessing()

        dfs = {}

        for sample in self.samples:
            print(f"Computing coverage: {sample.id}")
            df = sample.coverage_table()
            dfs[sample.id] = df

        z_df = self.compute_zscore(dfs)

        z_df.to_csv("zscore.tsv", sep="\t", index=False)


if __name__ == "__main__":
    pipeline = Pipeline(
        samples_table="samples.tsv",
        reference="reference.fa"
    )

    pipeline.run()