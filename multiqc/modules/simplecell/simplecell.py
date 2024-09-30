import logging
import csv
from typing import Dict, Union

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph

from .cell_caller import parse_cellcaller_json
from .cell_caller_mixed_species import parse_cellcaller_mixed_json

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module parses summary statistics from the `<sample_id>_metrics.csv` output files.
    Sample names are taken from the filename prefix.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="simplecell",
            anchor="simplecell",
            href="https://github.com/csgenetics/csgenetics_scrnaseq",
            info="CS Genetics simplecell scRNA-seq processing pipeline",
        )

        # Find all files for simplecell
        self.data_by_sample: Dict[str, int] = dict()
        for f in self.find_log_files("simplecell/metrics"):
            self.add_data_source(f)
            # TODO: log simplecell version
            s_name = f["s_name"]
            self.add_software_version(None, s_name)
            if s_name in self.data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.data_by_sample[s_name] = self.parse_file(f["f"])

        headers = {
            "num_cells": {
                "title": "Number of cells",
                "description": "Estimated number of cells; Number of barcodes passing the total counts threshold",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "raw_reads_per_cell": {
                "title": "Raw reads per cell",
                "description": "Number of reads pre-QC / Number of cells",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "median_genes_detected_per_cell": {
                "title": "Median genes detected per cell",
                "description": "Median number of genes detected for the cells (including nuclear and mitochondrial genes)",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        }
        # Select metrics for general stats table
        general_stats_metrics = ["num_cells", "raw_reads_per_cell", "median_genes_detected_per_cell"]
        self.general_stats_data = {
            sample: {key: value for key, value in metrics.items() if key in general_stats_metrics}
            for sample, metrics in self.data_by_sample.items()
        }
        # Add data to general stats table
        self.general_stats_addcols(self.general_stats_data, headers)
        # Write data to multiqc out file
        self.write_data_file(self.data_by_sample, "simplecell_metrics")
        log.info(f"Found {len(self.data_by_sample)} reports")

        # Add plot/table sections
        ## Define headers for read qc table
        ## Only metrics with defined header will be included in table
        read_qc_header = {
            "reads_pre_qc": {
                "title": "Number of reads pre-QC",
                "description": "Number of reads in the input R1 fastq files (after merging if applicable).",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "valid_barcode_reads": {
                "title": "Number of valid barcode-containing reads",
                "description": "Number of reads containing a barcode exactly matching the barcode_list.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "valid_barcode_reads_perc": {
                "title": "Percentage valid barcode-containing reads",
                "description": "(Number of valid barcode-containing reads / Number of reads pre-QC) * 100.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "barcode_bases_q30_perc": {
                "title": "Barcode bp >= Q30 percentage",
                "description": "The percentage of the barcode bases with a Phred score >= 30.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "reads_post_trimming": {
                "title": "Number of reads post-QC trimming",
                "description": "Number of reads after polyX tail and polyA internal trimming.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "reads_post_trimming_perc": {
                "title": "Percentage reads post-QC trimming",
                "description": "(Number of reads after polyX tail and polyA internal trimming / Number of valid barcode-containing reads) * 100.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "mean_post_trim_read_length": {
                "title": "Mean read length post-QC trimming",
                "description": "Mean R1 read length post-QC trimming.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "rna_bases_q30_perc": {
                "title": "R1 bp >= Q30 percentage; post-QC trimming",
                "description": "The percentage of the R1 bases (post-QC trimming) with a Phred score >= 30.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        }
        table_config = {"namespace": "", "id": "read-qc", "title": "Read QC", "only_defined_headers": True}
        # Add read QC section
        self.add_section(
            name="Read QC", plot=table.plot(self.data_by_sample, headers=read_qc_header, pconfig=table_config)
        )

        # Deduplication data
        for sample, metrics in self.data_by_sample.items():
            if "reads_before_deduplication" in metrics and "reads_after_deduplication" in metrics:
                metrics["reads_lost_deduplication"] = (
                    metrics["reads_before_deduplication"] - metrics["reads_after_deduplication"]
                )
        # Deduplication section headers
        cats = {
            "reads_after_deduplication": {
                "name": "Reads after deduplication",
                "description": "Number of high confidence (unique alignment: max mismatch <= 3bp) gene-annotated (with XT gene_id annotation) reads after deduplication.",
                "color": "#8bbc21",
            },
            "reads_lost_deduplication": {
                "name": "Reads lost during deduplication",
                "description": "Reads before deduplication - reads after deduplication",
                "color": "#f7a35c",
            },
        }
        dedup_config = {"id": "<Read deduplication>", "title": "Read deduplication"}
        self.add_section(name="Deduplication", plot=bargraph.plot(self.data_by_sample, cats, pconfig=dedup_config))

        # Add cell metrics section
        cell_metrics_header = {
            "median_total_reads_per_cell": {
                "title": "Median total counts per cell",
                "description": ",Median sum of counts per cell.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "median_genes_detected_per_cell": {
                "title": "Median genes detected per cell",
                "description": "Median sum of counts per cell.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "median_nuclear_genes_detected_per_cell": {
                "title": "Median nuclear (non-mitochondrial) genes detected per cell",
                "description": "Median number of nuclear genes detected for each cell.",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        }
        table_config = {"namespace": "", "id": "cell-metrics", "title": "Cell metrics", "only_defined_headers": True}
        # Add read QC section
        self.add_section(
            name="Cell metrics", plot=table.plot(self.data_by_sample, headers=cell_metrics_header, pconfig=table_config)
        )

        parse_cellcaller_json(self)
        parse_cellcaller_mixed_json(self)

    def parse_file(self, f) -> Dict[str, Union[float, int]]:
        data = {}
        reader = csv.reader(f.splitlines())
        next(reader)  # Skip the header row
        for row in reader:
            key = row[0].strip()
            value = float(row[1].strip())
            data[key] = value
        return data
