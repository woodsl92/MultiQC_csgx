import logging
from copy import deepcopy
import csv
from typing import Dict, Union, List

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots.table_object import ColumnMeta
from multiqc.plots import table

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
        self.data_by_sample: Dict[str, Dict[str, Union[float, int]]] = dict()
        for f in self.find_log_files("simplecell"):
            self.add_data_source(f)
            # TODO: log simplecell version
            s_name = f["s_name"]
            self.add_software_version(None, s_name)
            if s_name in self.data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.data_by_sample[s_name] = self.parse_file(f["f"])

        headers = {
            "num_cells": {
                "title": "Num cells",
                "description": "Estimated number of cells; Number of barcodes passing the total counts threshold",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "median_total_reads_per_cell": {
                "title": "Median reads per cell",
                "description": "Median sum of counts per cell",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
        }
        # Select metrics for general stats table
        general_stats_metrics = ["num_cells", "median_total_reads_per_cell"]
        self.general_stats_data = {
            sample: {key: value for key, value in metrics.items() if key in general_stats_metrics}
            for sample, metrics in self.data_by_sample.items()
        }
        # Add data to general stats table
        self.general_stats_addcols(self.general_stats_data, headers)
        # Write data to multiqc out file
        self.write_data_file(self.data_by_sample, "multiqc_simplecell")
        log.info(f"Found {len(self.data_by_sample)} reports")

        # Select metrics for csgenetics section table
        cs_stats_metrics = ["num_cells", "raw_reads_per_cell", "median_total_reads_per_cell"]
        self.cs_stats_data = {
            sample: {key: value for key, value in metrics.items() if key in cs_stats_metrics}
            for sample, metrics in self.data_by_sample.items()
        }
        single_header = {}
        table_config = {"namespace": "", "id": "per-cell-metrics", "title": "Per-cell metrics"}
        # Add cs table section
        self.add_section(plot=table.plot(self.cs_stats_data, headers=single_header, pconfig=table_config))

    def parse_file(self, f) -> Dict[str, Union[float, int]]:
        data = {}
        reader = csv.reader(f.splitlines())
        next(reader)  # Skip the header row
        for row in reader:
            key = row[0].strip()
            value = float(row[1].strip())
            data[key] = value
        return data
