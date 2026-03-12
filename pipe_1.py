class Sample:
    
    def __init__(self, sample_id, r1_path, r2_path):
        self.id = sample_id # 
        self.r1 = Path(r1_path) # путь до форварда
        self.r2 = Path(r2_path) # путь до обратного 
        self.sample_dir = output_dir / self.id # общая папка
        
        # папки, чтобы выкидывать туда результаты
        self.trimmed_dir = Path("trimmed")
        self.bam_dir = Path("bam")
        self.plots_dir = Path("plots")
        for d in [self.trimmed_dir, self.bam_dir, self.plots_dir]:
            d.mkdir(parents=True, exist_ok=True) # если есть, то кладет, если нет, то создает без ошибок
            
        # новые файлы, которые появятся в папках, инициализированных выше
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

    def coverage_table(self):

        depth_file = self.bam_dir / f"{self.id}.depth.tsv"

        with pysam.AlignmentFile(self.bam_file, "rb") as bam, open(depth_file, "w") as f:

            for col in bam.pileup():
                chrom = col.reference_name
                pos = col.pos + 1
                depth = col.nsegments
                line = f"{chrom}\t{pos}\t{depth}\n"
                f.write(line)

        df = pd.read_csv(depth_file, sep="\t", names=["chr", "pos", "depth"])

        return df
        
    def coverage_table_mb(self, window_size=1_000_000):
        depth_file = self.bam_dir / f"{self.id}.1mb.cov.tsv"

        with pysam.AlignmentFile(self.bam_file, "rb") as bam, open(depth_file, "w") as f:
            coverage = {}
            for col in bam.pileup():
                chrom = col.reference_name
                window = col.pos // window_size
                if (chrom, window) not in coverage:
                    coverage[(chrom, window)] = 0
                coverage[(chrom, window)] += col.nsegments

            for (chrom, window), depth_sum in sorted(coverage.items()):
                start = window * window_size
                end = start + window_size
                mean_depth = depth_sum / window_size
                f.write(f"{chrom}\t{start}\t{end}\t{mean_depth}\n")

    df = pd.read_csv(depth_file, sep="\t", names=["chr", "start", "end", "coverage"])
    return df
    
    def plot_chromosome_coverage(self, chrom="chr1", window_size=100_000):
        
        df = pd.read_csv(
            self.coverage_file,
            sep="\t",
            header=None,
            names=["chr", "pos", "depth"]
        )

        df_chr = df[df["chr"] == chrom].copy()

        df_chr["window"] = df_chr["pos"] // window_size

        windowed = df_chr.groupby("window")["depth"].mean().reset_index()
        windowed["pos"] = windowed["window"] * window_size

        plt.figure(figsize=(20, 5))
        plt.plot(windowed["pos"], windowed["depth"], color="blue")
        plt.xlabel("Position on " + chrom)
        plt.ylabel("Average coverage per 100kb")
        plt.title(f"Smoothed coverage on {chrom} ({self.sample_id})")
        plt.grid(True)
        plt.tight_layout()
        plt.show()


class Pipeline:
    def __init__(self, samples_table, reference):
        self.samples_table = Path(samples_table)
        self.reference = Path(reference)
        self.samples = self.load_samples()

    def load_samples(self):
        df = pd.read_csv(self.samples_table, sep="\t")
        samples = []

        for _, row in df.iterrows():
            sample = Sample(
                sample_id=row["ID"],
                r1_path=row["R1"],
                r2_path=row["R2"],
                output_dir=self.output_dir
            )
            samples.append(sample)
        return samples
        
    def compute_zscore(self, dfs, total_reads):
        tables = []
    
        for sample_id, df in dfs.items():
            df = df[["chr","start","coverage"]].copy()
            df[sample_id] = df["coverage"] / total_reads[sample_id] * 1_000_000
            tables.append(df[["chr","start",sample_id]])

        merged = tables[0]

        for x in tables[1:]:
            merged = merged.merge(t, on=["chr","start"])

        sample_cols = list(dfs())

        mean = merged[sample_cols].mean(axis=1)
        std = merged[sample_cols].std(axis=1)

        for col in sample_cols:
            merged[col+"_z"] = (merged[col] - mean) / std

        return merged

    def plot_chromosomes(self, df):
        chrom_means = df.groupby("chr").mean(numeric_only=True)

        plt.figure(figsize=(12,5))

        for col in df.columns:
            if col.endswith("_z"):
                plt.scatter(
                    chrom_means.index,
                    chrom_means[col],
                    s=80,
                    label=col.replace("_z",""))

        plt.xlabel("Chromosome")
        plt.ylabel("Mean Z-score")
        plt.title("Mean Z-score per chromosome")

        plt.xticks(rotation=45)
        plt.legend()
        plt.tight_layout()
        plt.show()
            
    def plot_bins(self, df):
        x = np.arange(len(df))

        plt.figure(figsize=(15,6))

        for col in df.columns:
            if col.endswith("_z"):
                plt.scatter(
                    x,
                    df[col],
                    s=3,
                    alpha=0.6,
                    label=col.replace("_z",""))

        plt.axhline(0, linestyle="--")
        plt.xlabel("Genomic bins")
        plt.ylabel("Z-score")
        plt.title("Z-score per genomic bin")

        plt.legend()
        plt.tight_layout()
        plt.show()  
            
    def run(self):
        for sample in self.samples:
            sample.trim()
            sample.align(self.reference)

            df = sample.coverage_table_mb()

            dfs[sample.id] = df
            total_reads[sample.id] = df["coverage"].sum()

        z_df = self.compute_zscore(dfs, total_reads)

        self.plot_chromosomes(z_df)
        self.plot_bins(z_df)
        return z_df