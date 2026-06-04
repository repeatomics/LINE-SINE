class MetricsCollector:

    """
    Container for collecting per-sample QC metrics.

    Stores metrics in memory and exports them as a TSV table.
    """

    def __init__(self):
        """Initialize empty metrics storage."""
        self.rows = []
        
    def add(self, sample_id, metric, value):
        """
        Add a metric value for a sample.

        Args:
            sample_id -> str: sample name
            metric -> str: metric name 
            value -> float: numeric value 
        """
        self.rows.append({
            "sample": sample_id,
            "metric": metric,
            "value": value
        })

    def save(self, path):
        """
        Save all metrics to TSV file.

        Args:
            path -> str: output path
        """
        pd.DataFrame(self.rows).to_csv(path, sep="\t", index=False)
