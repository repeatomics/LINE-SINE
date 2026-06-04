from pathlib import Path
from sample import Sample
from metrics import MetricsCollector
from config import Config


class Pipeline:
    """
    Main controller of the sequencing pipeline.

    Runs full analysis workflow from raw data to final metrics:

    Steps:
     Sample detection (FASTQ or BAM mode)
     Trimming (fastp)
     Alignment (BWA-MEM)
     Coverage calculation
     Variant calling (bcftools)
     VAF computation
     Metrics comparation

    Input:
        Directory with FASTQ or BAM files + reference genome

    Output:
     metrics.tsv 
     per-sample coverage data
     VCF files
     VAF tables
    """
    def __init__(self, root, ref, config=None):

        self.root = Path(root)
        self.ref = ref

        self.config = config or Config()

        self.samples = []
        self.metrics = MetricsCollector()

        self.find_samples()

    # FIND
    def find_samples(self):

        run_name = self.root.name

        r1_files = list(self.root.rglob("*R1*.fastq.gz"))

        if len(r1_files) > 0:

            for r1 in r1_files:

                r2 = Path(str(r1).replace("_R1", "_R2"))
                if not r2.exists():
                    continue

                sid = r1.name.split("_R1")[0]

                self.samples.append(
                    Sample(sid, run_name, r1=r1, r2=r2)
                )

        else:

            bam_files = list(self.root.rglob("*.bam"))

            for bam in bam_files:

                sid = bam.stem.replace(".sorted", "")

                self.samples.append(
                    Sample(sid, run_name, bam=bam)
                )

    # RUN
    def run(self):

        for s in self.samples:

            print(f"\n[{s.id}] START")

            if not Path(s.bam).exists():

                s.trim(self.config)
                s.align(self.ref, self.config)

            print(f"[{s.id}] coverage")

            cov = s.coverage(
                self.config.bin_size,
                self.config.min_mapping_quality
            )

            if cov.empty:
                continue

            self.metrics.add(s.id, "mean_cov", cov["coverage"].mean())

            s.plot_genome_coverage(cov)

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