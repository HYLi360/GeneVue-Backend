# coding: utf-8
"""
BuscoPlot.py

Introduced in v6.0.0, the BuscoPlot module is a rewrite of the previous python/R script.
The plotting functionality is fully integrated into the BUSCO code base and is accessed
via the `BUSCO --plot` command.

Author(s): Matthew Berkeley

Copyright (c) 2015-2025, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

import matplotlib.pyplot as plt
import numpy as np
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
import os
import json

logger = BuscoLogger.get_logger(__name__)


class Plot:
    def __init__(self, config, wd=None):
        if wd:
            self.results_dir = wd
        else:
            if "/" in config.get("busco_run", "out"):
                self.results_dir = os.path.expanduser(config.get("busco_run", "out"))
            else:
                self.results_dir = os.path.join(
                    config.get("busco_run", "out_path"), config.get("busco_run", "out")
                )

        self.check_wd(self.results_dir)

        self.data = {
            "labels": [],
            "values": [],
            "percentages": [],
            "datasets": [],
            "totals": [],
            "summaries": [],
        }

        self.categories = [
            "Complete (C) and single-copy (S)",
            "Complete (C) and duplicated (D)",
            "Fragmented (F)",
            "Missing (M)",
        ]
        self.colors = ["#56B4E9", "#3492C7", "#F0E442", "#F04442"]

        self.plot_percentages = config.getboolean("busco_run", "plot_percentages")

    @staticmethod
    def check_wd(wd):
        """
        This function checks that the working directory exists with write permission
        :raises SystemExit: if the folder is absent or the user has no write permission
        """
        if not os.path.exists(wd):
            logger.warning("Directory {} does not exist".format(wd))
            raise SystemExit()
        if not os.access(wd, os.W_OK):
            logger.warning("No permission to write into {}".format(wd))
            raise SystemExit()

    def generate_plot(self):
        max_labels_per_page = 10
        num_pages = (
            len(self.data["labels"]) + max_labels_per_page - 1
        ) // max_labels_per_page

        for page in range(num_pages):
            # Determine the range of labels for the current page
            start_idx = page * max_labels_per_page
            end_idx = min(start_idx + max_labels_per_page, len(self.data["labels"]))

            # Extract data for the current page
            labels = self.data["labels"][start_idx:end_idx]
            percentages = self.data["percentages"][start_idx:end_idx]
            summaries = self.data["summaries"][start_idx:end_idx]

            # Plot
            fig, ax = plt.subplots(figsize=(10, 6))
            bar_width = 0.5
            y_positions = np.arange(len(labels))

            percentages_array = np.array(percentages)  # Convert to NumPy array
            # Stacked bar plot
            for i, category in enumerate(self.categories):
                bottom = percentages_array[:, 1 : i + 1].sum(axis=1) if i > 0 else None
                ax.barh(
                    y_positions,
                    percentages_array[:, i + 1],
                    color=self.colors[i],
                    label=category,
                    left=bottom,
                    height=bar_width,
                )

            # Add one_line_summary text on top of each bar group
            for y_pos, summary in zip(y_positions, summaries):
                fontsize = max(
                    6, min(12, bar_width * 10 / len(labels))
                )  # Scale fontsize
                ax.text(
                    5,  # Position at the left end of the stacked bar
                    y_pos,
                    summary,
                    ha="left",
                    va="center",
                    fontsize=fontsize,
                    color="black",
                    fontweight="bold",
                )

            # Customization
            ax.set_yticks(y_positions)
            ax.set_yticklabels(labels)
            ax.set_xlabel("% BUSCOs")
            ax.set_title(f"BUSCO Assessment Results (Page {page + 1})")
            ax.legend(
                loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False
            )
            ax.grid(axis="x", linestyle="--", alpha=0.7)

            # Add space between the y-axis and the bar zero point
            min_value = percentages_array.min() if percentages_array.min() < 0 else 0
            ax.set_xlim(left=min_value - 5)  # Add a margin of 5 units to the left

            # Adjust layout and save
            plt.tight_layout()
            output_filename = (
                f"busco_figure_page_{page + 1}.png"
                if num_pages > 1
                else "busco_figure.png"
            )
            plt.savefig(
                os.path.join(self.results_dir, output_filename),
                dpi=300,
            )
            plt.close(fig)

    def load_single_json(self, json_file):

        with open(json_file, "r") as content:
            results_data = json.load(content)
        dataset = results_data["lineage_dataset"]["name"]
        dataset_creation_date = results_data["lineage_dataset"]["creation_date"]
        busco_version = results_data["versions"]["busco"]
        input_file = os.path.basename(
            results_data["parameters"]["in"]
        )  # remove full path
        comp = results_data["results"]["Complete BUSCOs"]
        comp_pc = results_data["results"]["Complete percentage"]
        sc = results_data["results"]["Single copy BUSCOs"]
        sc_pc = results_data["results"]["Single copy percentage"]
        dupl = results_data["results"]["Multi copy BUSCOs"]
        dupl_pc = results_data["results"]["Multi copy percentage"]
        frag = results_data["results"]["Fragmented BUSCOs"]
        frag_pc = results_data["results"]["Fragmented percentage"]
        miss = results_data["results"]["Missing BUSCOs"]
        miss_pc = results_data["results"]["Missing percentage"]
        total = results_data["results"]["n_markers"]
        one_line_summary = "C:{}[S:{},D:{}],F:{},M:{},n:{}".format(
            comp, sc, dupl, frag, miss, total
        )
        one_line_summary_pc = results_data["results"]["one_line_summary"]

        if "plot_label" in results_data["parameters"]:
            self.data["labels"] += [results_data["parameters"]["plot_label"]]
        else:
            self.data["labels"] += [input_file]

        if self.plot_percentages:
            self.data["values"].append(
                [
                    comp_pc,
                    sc_pc,
                    dupl_pc,
                    frag_pc,
                    miss_pc,
                ]
            )
            summary_line = one_line_summary_pc

        else:
            self.data["values"].append([comp, sc, dupl, frag, miss])
            summary_line = one_line_summary

        self.data["summaries"] += [summary_line]

        self.data["percentages"].append(
            [
                comp_pc,
                sc_pc,
                dupl_pc,
                frag_pc,
                miss_pc,
            ]
        )

        self.data["datasets"] += [dataset]
        self.data["totals"] += [total]

        logger.info("Loaded {} successfully".format(json_file))

    def load_data(self):

        if "batch_summary.txt" in os.listdir(self.results_dir):
            logger.info("Plotting batch results from {}".format(self.results_dir))
            results_dirs = [
                os.path.join(self.results_dir, d)
                for d in os.listdir(self.results_dir)
                if os.path.isdir(os.path.join(self.results_dir, d)) and d != "logs"
            ]
            standard_json_filenames = []
            for resdir in results_dirs:
                ## file that begins with "short_summary" and ends with ".json"
                standard_json_filenames += [
                    os.path.join(resdir, f)
                    for f in os.listdir(resdir)
                    if f.startswith("short_summary") and f.endswith(".json")
                ]
        else:
            # load data from json files
            standard_json_filenames = [
                os.path.join(self.results_dir, f)
                for f in os.listdir(self.results_dir)
                if f.endswith(".json")
            ]
        if len(standard_json_filenames) > 0:
            for f in standard_json_filenames:
                try:
                    self.load_single_json(f)
                except IOError:
                    logger.warning("Impossible to use the file {}".format(f))
                except json.decoder.JSONDecodeError:
                    logger.warning("Impossible to read the JSON file {}".format(f))
        else:
            logger.warning(
                "No short_summary JSON files found in {}".format(self.results_dir)
            )

        return
