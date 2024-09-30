"""Parse output from simplecell cellcaller mixed species"""

import json
import logging
from typing import Dict
from multiqc import BaseMultiqcModule, config
from multiqc.plots import linegraph


log = logging.getLogger(__name__)


def parse_cellcaller_mixed_json(module: BaseMultiqcModule) -> int:
    """
    Simplecell cellcaller plot data parser
    """
    # Define the nested dictionary type
    # Initialize the data_by_sample dictionary with the correct type annotation
    data_by_sample: Dict[str, Dict[str, Dict]] = dict()
    hsap_data_by_sample: Dict[str, Dict] = dict()
    hsap_noise_by_sample = {}
    hsap_cells_by_sample = {}
    hsap_noise_cells_by_sample = {}

    mmus_noise_by_sample = {}
    mmus_cells_by_sample = {}
    mmus_noise_cells_by_sample = {}

    for f in module.find_log_files("simplecell/cellcaller_mixed_species", filehandles=True):
        module.add_data_source(f)
        s_name = f["s_name"]
        module.add_software_version(None, s_name)
        data_by_sample[s_name] = json.load(f["f"])
        # Extract hsap/mmus data

        # Extract noisy/cell barcode plotting data, update names to plot 2 lines per sample
        hsap_data_by_sample = data_by_sample[s_name]["hsap_pd_data"]
        hsap_noise_by_sample[f"{s_name}_noise"] = hsap_data_by_sample["noise"]
        hsap_cells_by_sample[f"{s_name}_cell"] = hsap_data_by_sample["cell"]

        mmus_data_by_sample = data_by_sample[s_name]["mmus_pd_data"]
        mmus_noise_by_sample[f"{s_name}_noise"] = mmus_data_by_sample["noise"]
        mmus_cells_by_sample[f"{s_name}_cell"] = mmus_data_by_sample["cell"]

    # Combine noise/cell dicts into noise_cell dict
    hsap_noise_cells_by_sample.update(hsap_noise_by_sample)
    hsap_noise_cells_by_sample.update(hsap_cells_by_sample)
    mmus_noise_cells_by_sample.update(mmus_noise_by_sample)
    mmus_noise_cells_by_sample.update(mmus_cells_by_sample)

    module.write_data_file(data_by_sample, "simplecell_cellcaller_mixed")
    log.info(f"Found {len(data_by_sample)} reports")

    ## Hsap plot
    # Get max value of y to plot default threshold line
    hsap_max_y = max(hsap_noise_cells_by_sample[f"{s_name}_noise"].values())
    hsap_max_cell = max(hsap_noise_cells_by_sample[f"{s_name}_cell"].values())
    if hsap_max_cell > hsap_max_y:
        hsap_max_y = hsap_max_cell

    hsap_line_config = {
        "id": "cellcaller_mixed_hsap",
        "title": "Cellcaller mixed hsap",
        "anchor": "cellcaller_mixed_hsap",
        "xlab": "log10(total counts+1)",
        "ylab": "Density",
        "showlegend": True,
        "extra_series": [
            {
                "name": "Default threshold",
                "pairs": [[2, 0], [2, hsap_max_y]],
                "dash": "dash",
                "width": 1,
                "color": "#000000",
            }
        ],
    }

    ## mmus plot
    # Get max value of y to plot default threshold line
    mmus_max_y = max(mmus_noise_cells_by_sample[f"{s_name}_noise"].values())
    mmus_max_cell = max(mmus_noise_cells_by_sample[f"{s_name}_cell"].values())
    if mmus_max_cell > mmus_max_y:
        mmus_max_y = mmus_max_cell

    mmus_line_config = {
        "id": "cellcaller_mixed_mmus",
        "title": "Cellcaller mixed mmus",
        "anchor": "cellcaller_mixed_mmus",
        "xlab": "log10(total counts+1)",
        "ylab": "Density",
        "showlegend": True,
        "extra_series": [
            {
                "name": "Default threshold",
                "pairs": [[2, 0], [2, mmus_max_y]],
                "dash": "dash",
                "width": 1,
                "color": "#000000",
            }
        ],
    }

    module.add_section(
        name="Cellcaller hsap",
        anchor="cellcaller-hsap",
        description="Plots produced by simplecell cellcaller - hsap",
        plot=linegraph.plot(hsap_noise_cells_by_sample, hsap_line_config),
    )

    module.add_section(
        name="Cellcaller mmus",
        anchor="cellcaller-mmus",
        description="Plots produced by simplecell cellcaller - mmus",
        plot=linegraph.plot(mmus_noise_cells_by_sample, mmus_line_config),
    )
