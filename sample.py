import subprocess
from pathlib import Path
import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt

from config import Config


class Sample:
    """
    A single sequencing sample in the pipeline.

    This object stores all sample-related files and runs processing steps:
    trimming, alignment, coverage calculation, variant calling, and VAF estimation.

    Input:
        FASTQ (R1/R2) or pre-aligned BAM
    Output:
        Processed BAM, coverage tables, VCF, VAF table, QC reports
    """

    def __init__(self, sample_id, run, r1=None, r2=None, bam=None):
        
        """
        Initialize a Sample.

        Args:
            sample_id -> str: sample name 
            run -> str: dataset or run name
            r1 -> path or str: FASTQ R1 file
            r2 --> path or str: FASTQ R2 file
            bam -> path or str: pre-aligned BAM file
        """


        self.id = str(sample_id)
        self.run = run

        self.r1 = r1
        self.r2 = r2
        self.bam_input = bam

        self.base = Path("output") / run / self.id

        self.trim_dir = self.base / "trim"
        self.bam_dir = self.base / "bam"
        self.plots_dir = self.base / "plots"
        self.qc_dir = self.base / "qc"

        for d in [self.trim_dir, self.bam_dir, self.plots_dir, self.qc_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.r1_trim = self.trim_dir / "R1.fastq.gz"
        self.r2_trim = self.trim_dir / "R2.fastq.gz"

        self.bam = self.bam_input or self.bam_dir / f"{self.id}.bam"

        self.coverage_file = self.bam_dir / f"{self.id}.cov.tsv"
        self.vaf_file = self.base / "vaf.tsv"

        self.fastp_html = self.qc_dir / "fastp.html"
        self.fastp_json = self.qc_dir / "fastp.json"

        if (
            self.bam_input
            and Path(self.bam).exists()
            and not Path(str(self.bam) + ".bai").exists()
        ):
            subprocess.run(["samtools", "index", str(self.bam)], check=True)

    # TRIM
    def trim(self, config: Config):

        if self.r1 is None or self.r2 is None:
            raise ValueError("FASTQ not provided")

        cmd = [
            "fastp",
            "-i", str(self.r1),
            "-I", str(self.r2),
            "-o", str(self.r1_trim),
            "-O", str(self.r2_trim),
            "--qualified_quality_phred", str(config.fastp_quality_phred),
            "--length_required", str(config.fastp_length_required),
            "-h", str(self.fastp_html),
            "-j", str(self.fastp_json),
        ]

        subprocess.run(cmd, check=True)

    # ALIGN
    def align(self, reference, config: Config):

        cmd = (
            f"bwa mem -t {config.bwa_threads} {reference} "
            f"{self.r1_trim} {self.r2_trim} | "
            f"samtools view -bS - | "
            f"samtools sort -@ {config.samtools_threads} -m 1G -o {self.bam}"
        )

        subprocess.run(cmd, shell=True, check=True)

        subprocess.run(["samtools", "index", str(self.bam)], check=True)

    # COVERAGE
    def coverage(self, bin_size, min_mapq):

        if not Path(self.bam).exists():
            return pd.DataFrame(columns=["chr", "start", "coverage"])

        cov = {}

        with pysam.AlignmentFile(self.bam, "rb") as bam:
            for col in bam.pileup(min_mapping_quality=min_mapq):

                chrom = col.reference_name
                if chrom is None:
                    continue

                window = col.reference_pos // bin_size
                key = (chrom, window)

                cov[key] = cov.get(key, 0) + col.nsegments

        rows = []

        with open(self.coverage_file, "w") as f:
            for (chrom, window), value in sorted(cov.items()):

                start = window * bin_size
                mean_cov = value / bin_size

                rows.append([chrom, start, mean_cov])
                f.write(f"{chrom}\t{start}\t{mean_cov}\n")

        return pd.DataFrame(rows, columns=["chr", "start", "coverage"])

    # PLOT
    def plot_genome_coverage(self, cov_df):

        if cov_df.empty:
            return

        x = np.arange(len(cov_df))

        plt.figure(figsize=(20, 5))
        plt.plot(x, cov_df["coverage"].values, linewidth=1)

        plt.xlabel("Genome bins")
        plt.ylabel("Coverage")
        plt.title(self.id)
        plt.grid(alpha=0.3)
        plt.tight_layout()

        out = self.plots_dir / "genome_coverage.png"
        plt.savefig(out, dpi=300)
        plt.close()

    # VARIANTS
    def call_variants(self, reference):

        vcf = self.base / f"{self.id}.vcf.gz"

        cmd = (
            f"bcftools mpileup -f {reference} -a AD {self.bam} | "
            f"bcftools call -mv -Oz -o {vcf}"
        )

        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(["bcftools", "index", str(vcf)], check=True)

        return vcf

    # VAF
    def compute_vaf(self, vcf):

        cmd = (
            f"bcftools query -f '%CHROM\t%POS\t[%AD]\n' {vcf}"
        )

        try:
            out = subprocess.check_output(cmd, shell=True).decode()
        except:
            return pd.DataFrame(columns=["chr", "pos", "vaf"])

        rows = []

        for line in out.strip().split("\n"):

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            chrom, pos, ad = parts

            if ad == "." or "," not in ad:
                continue

            try:
                ref, alt = ad.split(",")[:2]
                ref = float(ref)
                alt = float(alt)
            except:
                continue

            if ref + alt == 0:
                continue

            vaf = alt / (ref + alt)

            rows.append([chrom, int(pos), vaf])

        df = pd.DataFrame(rows, columns=["chr", "pos", "vaf"])
        df.to_csv(self.vaf_file, sep="\t", index=False)

        return df