"""Parse output from simplecell cellcaller single species"""

import json
import logging
from typing import Dict
from multiqc import BaseMultiqcModule, config
from multiqc.plots import linegraph


log = logging.getLogger(__name__)


def parse_cellcaller_json(module: BaseMultiqcModule) -> int:
    """
    Simplecell cellcaller plot data parser
    """
    data_by_sample: Dict[str, Dict] = dict()
    noise_by_sample = {}
    cells_by_sample = {}
    noise_cells_by_sample = {}

    for f in module.find_log_files("simplecell/cellcaller", filehandles=True):
        module.add_data_source(f)
        s_name = f["s_name"]
        module.add_software_version(None, s_name)
        data_by_sample[s_name] = json.load(f["f"])

        # Extract noisy/cell barcode plotting data, update names to plot 2 lines per sample
        noise_by_sample[f"{s_name}_noise"] = data_by_sample[s_name]["noise"]
        cells_by_sample[f"{s_name}_cell"] = data_by_sample[s_name]["cell"]
    # Combine noise/cell dicts into noise_cell dict
    noise_cells_by_sample.update(noise_by_sample)
    noise_cells_by_sample.update(cells_by_sample)

    module.write_data_file(data_by_sample, "simplecell_cellcaller")
    log.info(f"Found {len(data_by_sample)} reports")

    # Get max value of y to plot default threshold line
    max_y = max(noise_cells_by_sample[f"{s_name}_noise"].values())
    max_cell = max(noise_cells_by_sample[f"{s_name}_cell"].values())
    if max_cell > max_y:
        max_y = max_cell

    line_config = {
        "id": "cellcaller",
        "title": "Cellcaller",
        "anchor": "cellcaller",
        "xlab": "log10(total counts+1)",
        "ylab": "Density",
        "showlegend": True,
        "extra_series": [
            {
                "name": "Default threshold",
                "pairs": [[2, 0], [2, max_y]],
                "dash": "dash",
                "width": 1,
                "color": "#000000",
            }
        ],
    }
    module.add_section(
        name="Cellcaller",
        anchor="cellcaller",
        description="Plots produced by simplecell cellcaller",
        plot=linegraph.plot(noise_cells_by_sample, line_config),
    )