import pandas as pd
import pysam
from pathlib import Path
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import sys


# METRICS

class MetricsCollector:

    def __init__(self):
        self.rows = []

    def add(self, sample_id, metric, value):
        self.rows.append({
            "sample": sample_id,
            "metric": metric,
            "value": value
        })

    def save(self, path):
        pd.DataFrame(self.rows).to_csv(path, sep="\t", index=False)


# SAMPLE

class Sample:

    def __init__(self, sample_id, run, r1=None, r2=None, bam=None):

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

        for d in [
            self.trim_dir,
            self.bam_dir,
            self.plots_dir,
            self.qc_dir
        ]:
            d.mkdir(parents=True, exist_ok=True)

        # trimmed fastq
        self.r1_trim = self.trim_dir / "R1.fastq.gz"
        self.r2_trim = self.trim_dir / "R2.fastq.gz"

        # bam
        self.bam = self.bam_input or self.bam_dir / f"{self.id}.bam"

        # coverage
        self.coverage_file = self.bam_dir / f"{self.id}.1mb.cov.tsv"

        # VAF
        self.vaf_file = self.base / "vaf.tsv"

        # fastp reports
        self.fastp_html = self.qc_dir / "fastp.html"
        self.fastp_json = self.qc_dir / "fastp.json"

        # auto-index BAM
        if (
            self.bam_input
            and Path(self.bam).exists()
            and not Path(str(self.bam) + ".bai").exists()
        ):
            subprocess.run(
                ["samtools", "index", str(self.bam)],
                check=True
            )

    # TRIM

    def trim(self):

        if self.r1 is None or self.r2 is None:
            raise ValueError("FASTQ not provided")

        cmd = [
            "fastp",
            "-i", str(self.r1),
            "-I", str(self.r2),
            "-o", str(self.r1_trim),
            "-O", str(self.r2_trim),
            "--qualified_quality_phred", "20",
            "--length_required", "30",
            "-h", str(self.fastp_html),
            "-j", str(self.fastp_json)
        ]

        subprocess.run(cmd, check=True)

    # ALIGN

    def align(self, reference):

        cmd = (
            f"bwa mem -t 6 {reference} "
            f"{self.r1_trim} {self.r2_trim} | "
            f"samtools view -bS - | "
            f"samtools sort -@ 4 -m 1G -o {self.bam}"
        )

        subprocess.run(cmd, shell=True, check=True)

        subprocess.run(
            ["samtools", "index", str(self.bam)],
            check=True
        )

    # COVERAGE

    def coverage(self, bin_size=1_000_000):

        if not Path(self.bam).exists():
            return pd.DataFrame(
                columns=["chr", "start", "coverage"]
            )

        cov = {}

        with pysam.AlignmentFile(self.bam, "rb") as bam:

            for col in bam.pileup(min_mapping_quality=20):

                chrom = col.reference_name

                if chrom is None:
                    continue

                window = col.reference_pos // bin_size

                key = (chrom, window)

                if key not in cov:
                    cov[key] = 0

                cov[key] += col.nsegments

        rows = []

        with open(self.coverage_file, "w") as f:

            for (chrom, window), value in sorted(cov.items()):

                start = window * bin_size
                mean_cov = value / bin_size

                rows.append([
                    chrom,
                    start,
                    mean_cov
                ])

                f.write(
                    f"{chrom}\t{start}\t{mean_cov}\n"
                )

        return pd.DataFrame(
            rows,
            columns=["chr", "start", "coverage"]
        )

    # PLOT COVERAGE

    def plot_genome_coverage(self, cov_df):

        if cov_df.empty:
            return

        x = np.arange(len(cov_df))

        plt.figure(figsize=(20, 5))

        plt.plot(
            x,
            cov_df["coverage"].values,
            linewidth=1
        )

        plt.xlabel("Genome bins")
        plt.ylabel("Coverage")

        plt.title(f"{self.id}")

        plt.grid(alpha=0.3)

        plt.tight_layout()

        out = self.plots_dir / "genome_coverage.png"

        plt.savefig(out, dpi=300)

        plt.close()

    # VARIANT CALLING

    def call_variants(self, reference):

        vcf = self.base / f"{self.id}.vcf.gz"

        cmd = (
            f"bcftools mpileup "
            f"-f {reference} "
            f"-a AD "
            f"{self.bam} | "
            f"bcftools call -mv -Oz -o {vcf}"
        )

        subprocess.run(
            cmd,
            shell=True,
            check=True
        )

        subprocess.run(
            ["bcftools", "index", str(vcf)],
            check=True
        )

        return vcf

    # VAF

    def compute_vaf(self, vcf):

        cmd = (
            f"bcftools query "
            f"-f '%CHROM\t%POS\t[%AD]\n' "
            f"{vcf}"
        )

        try:
            out = subprocess.check_output(
                cmd,
                shell=True
            ).decode()

        except:
            return pd.DataFrame(
                columns=["chr", "pos", "vaf"]
            )

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

            rows.append([
                chrom,
                int(pos),
                vaf
            ])

        df = pd.DataFrame(
            rows,
            columns=["chr", "pos", "vaf"]
        )

        df.to_csv(
            self.vaf_file,
            sep="\t",
            index=False
        )

        return df


# PIPELINE

class Pipeline:

    def __init__(self, root, ref):

        self.root = Path(root)
        self.ref = ref

        self.samples = []

        self.metrics = MetricsCollector()

        self.find_samples()

    # FIND SAMPLES

    def find_samples(self):

        run_name = self.root.name

        # FASTQ mode
        r1_files = list(
            self.root.rglob("*R1*.fastq.gz")
        )

        if len(r1_files) > 0:

            print("FASTQ MODE")
            print("FOUND:", len(r1_files))

            for r1 in r1_files:

                r2 = Path(
                    str(r1).replace("_R1", "_R2")
                )

                if not r2.exists():
                    continue

                sid = r1.name.split("_R1")[0]

                self.samples.append(
                    Sample(
                        sid,
                        run_name,
                        r1=r1,
                        r2=r2
                    )
                )

        # BAM mode
        else:

            bam_files = list(
                self.root.rglob("*.bam")
            )

            print("BAM MODE")
            print("FOUND:", len(bam_files))

            for bam in bam_files:

                sid = bam.stem.replace(".sorted", "")

                self.samples.append(
                    Sample(
                        sid,
                        run_name,
                        bam=bam
                    )
                )

    # RUN

    def run(self):

        print("SAMPLES:", len(self.samples))

        for s in self.samples:

            print(f"\n[{s.id}] START")

            # FASTQ -> ALIGN
            if not Path(s.bam).exists():

                print(f"[{s.id}] trim")
                s.trim()

                print(f"[{s.id}] align")
                s.align(self.ref)

            # COVERAGE
            print(f"[{s.id}] coverage")

            cov = s.coverage()

            if cov.empty:
                print(f"[{s.id}] no coverage")
                continue

            self.metrics.add(
                s.id,
                "mean_cov",
                cov["coverage"].mean()
            )

            s.plot_genome_coverage(cov)

            # VARIANTS
            print(f"[{s.id}] variants")

            vcf = s.call_variants(self.ref)

            vaf_df = s.compute_vaf(vcf)

            if not vaf_df.empty:

                self.metrics.add(
                    s.id,
                    "mean_vaf",
                    vaf_df["vaf"].mean()
                )

                self.metrics.add(
                    s.id,
                    "n_variants",
                    len(vaf_df)
                )

            print(f"[{s.id}] done")

        self.metrics.save("metrics.tsv")

        print("\nDONE")


# ENTRY POINT

if __name__ == "__main__":

    root = sys.argv[1]
    ref = sys.argv[2]

    Pipeline(root, ref).run()