from dataclasses import dataclass

@dataclass
class Config:
    """
    Pipeline configuration.

    Stores all user parameters.
    """

    bwa_threads: int = 6

    samtools_threads: int = 4

    min_mapping_quality: int = 20

    bin_size: int = 1_000_000

    fastp_quality_phred: int = 20

    fastp_length_required: int = 30