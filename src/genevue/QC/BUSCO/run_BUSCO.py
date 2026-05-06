#!/usr/bin/env python3
# coding: utf-8
"""
run_BUSCO.py

Main run script.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2025, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

import argparse
from argparse import RawTextHelpFormatter
import genevue.QC.BUSCO as BUSCO
from genevue.QC.BUSCO.BuscoRunner import AnalysisRunner, BatchRunner, SingleRunner
from genevue.QC.BUSCO.Exceptions import BatchFatalError, BuscoError, BuscoConfigError
from genevue.QC.BUSCO.base import ToolException
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
from genevue.QC.BUSCO.BuscoLogger import LogDecorator as log
from genevue.QC.BUSCO.ConfigManager import BuscoConfigManager
from genevue.QC.BUSCO.Actions import (
    ListLineagesAction,
    CleanHelpAction,
    CleanVersionAction,
    DirectDownload,
    GeneratePlotAction,
)
from genevue.QC.BUSCO.ConfigManager import BuscoConfigMain

import sys
import os
import time
import traceback
import socket
import requests
import re
import hashlib

logger = BuscoLogger.get_logger(__name__)


@log(
    "Start a [blue]BUSCO v{}[/blue] analysis, current time: [green]{}[/green]".format(
        BUSCO.__version__, time.strftime("%m/%d/%Y %H:%M:%S")
    ),
    logger,
)
class BuscoMaster:

    run_stats = {"versions": {"python": sys.version_info, "BUSCO": BUSCO.__version__}}

    def __init__(self, params):
        self.params = params
        self.config_manager = BuscoConfigManager(self.params)
        self.config = self.config_manager.config_main

    def load_config(self):
        """
        Load a BUSCO config file that will figure out all the params from all sources
        i.e., provided config file, dataset cfg, and user args
        """
        try:
            self.config_manager.load_busco_config_main()
        finally:
            self.config = self.config_manager.config_main

    def check_batch_mode(self):
        return self.config.getboolean("busco_run", "batch_mode")

    def run(self):
        try:
            self.load_config()
            BuscoLogger.init_file_handler()
            type(self).run_stats["config"] = self.config_manager.run_stats

            if self.config.getboolean(
                "busco_run", "restart"
            ):  # backup log file from last run
                try:
                    AnalysisRunner.backup_log_file(self.config)
                except:
                    pass
            runner = (
                BatchRunner(self.config_manager)
                if self.check_batch_mode()
                else SingleRunner(self.config_manager)
            )
            runner.run()
            type(self).run_stats.update(runner.run_stats)

        except BuscoConfigError as e:
            logger.error("Error in BUSCO configuration: {}".format(e))
            type(self).run_stats["config"] = self.config_manager.run_stats
            anon_str = re.sub(r"(/[^ \n\t]+)+", "/anonymised/filepath", str(e))
            type(self).run_stats["config"].update({"config_error": anon_str})
            raise SystemExit(1)

        except BuscoError as be:
            self.log_error(be)
            type(self).run_stats["error"] = str(be)
            raise SystemExit(1)

        except BatchFatalError as bfe:
            self.log_error(bfe)
            type(self).run_stats["error"] = str(bfe)
            raise SystemExit(1)

        except ToolException as e:
            self.log_error(e)
            type(self).run_stats["error"] = str(e)
            raise SystemExit(1)

        except KeyboardInterrupt:
            type(self).run_stats["error"] = "KeyboardInterrupt"
            self.log_error(
                "A signal was sent to kill the process. \nBUSCO analysis failed !"
            )
            raise SystemExit(1)

        except BaseException:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            logger.critical(
                "Unhandled exception occurred:\n{}\n".format(
                    "".join(
                        traceback.format_exception(exc_type, exc_value, exc_traceback)
                    )
                )
            )
            type(self).run_stats["error"] = str(exc_value)
            self.log_error(str(exc_value))
            raise SystemExit(1)

        finally:
            try:
                if not self.config.getboolean("busco_run", "opt-out-run-stats"):
                    logger.info(
                        "Thank you for using BUSCO! Anonymous usage data is gathered to improve the tool. You may opt out with --opt-out-run-stats."
                    )
                    self.assemble_run_data()
                    self.post_run_data()
            except:
                pass

            try:
                AnalysisRunner.move_log_file(self.config)
            except:
                pass
            finally:
                BuscoLogger.close_file_handler()

    @staticmethod
    def log_error(err):
        logger.error(err)
        logger.debug(err, exc_info=True)
        logger.error("BUSCO analysis failed!")
        logger.error(
            "Check the logs, read the user guide (https://busco.ezlab.org/busco_userguide.html), "
            "and check the BUSCO issue board on https://gitlab.com/ezlab/busco/issues\n"
        )

    def assemble_run_data(self):
        try:
            self.get_download_url()
        except:
            pass

        try:
            self.get_dist_info()
        except:
            pass

        try:
            self.get_ip_hash()
        except:
            pass

    def post_run_data(self):
        datetime_randomno = (
            time.strftime("%Y%m%d_%H%M%S_") + str(time.time()).split(".")[1]
        )

        headers = {"Content-Type": "application/json"}
        url = "https://busco-data.ezlab.org/uploaddb"
        type(self).run_stats["run_id"] = "".join(datetime_randomno.split("_"))

        response = requests.post(url, json=type(self).run_stats, headers=headers)

        if 200 <= response.status_code < 300:
            logger.debug("Data uploaded successfully.")
        else:
            logger.debug("Data upload failed. Status code: {}".format(response.status))

    def get_dist_info(self):
        if os.path.exists("/.dockerenv"):
            type(self).run_stats["distribution"] = "docker"
        elif os.path.exists("/.singularity.d"):
            type(self).run_stats["distribution"] = "singularity"
        else:
            type(self).run_stats["distribution"] = "manual"

    def get_download_url(self):
        if self.config.getboolean("busco_run", "offline"):
            type(self).run_stats["download_url"] = "offline"
        else:
            type(self).run_stats[
                "download_url"
            ] = self.config.downloader.download_base_url

    def get_ip_hash(self):
        try:
            hostname = socket.gethostname()
            ip_address = socket.gethostbyname(hostname)
            try:
                response = requests.get(f"https://ipinfo.io/{ip_address}/json")
                if response.status_code == 200:
                    data = response.json()
                    type(self).run_stats["country"] = data.get("country", "unknown")
            except:
                pass
            ip_hash = hashlib.md5(
                (ip_address + time.strftime("%m%Y")).encode()
            ).hexdigest()  # ensure anonymity by changing hash monthly
            type(self).run_stats["ip_hash"] = ip_hash
        except:
            pass


@log("Command line: {}".format(" ".join(sys.argv[:])), logger, debug=True)
def _parse_args():
    """
    This function parses the arguments provided by the user
    :return: a dictionary having a key for each argument
    :rtype: dict
    """

    parser = argparse.ArgumentParser(
        description="Welcome to BUSCO {}: the Benchmarking Universal Single-Copy Ortholog assessment tool.\n"
        "For more detailed usage information, please review the README file provided with "
        "this distribution and the BUSCO user guide. "
        "Visit this page https://gitlab.com/ezlab/busco#how-to-cite-busco to see how to cite BUSCO".format(
            BUSCO.__version__
        ),
        usage="BUSCO -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]",
        formatter_class=RawTextHelpFormatter,
        add_help=False,
    )

    optional = parser.add_argument_group("optional arguments")

    optional.add_argument(
        "-i",
        "--in",
        dest="in",
        required=False,
        metavar="SEQUENCE_FILE",
        help="Input sequence file in FASTA format. "
        "Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. "
        "Also possible to use a path to a directory containing multiple input files.",
    )

    optional.add_argument(
        "-o",
        "--out",
        dest="out",
        required=False,
        metavar="OUTPUT",
        help="Give your analysis run a recognisable short name. "
        "Output folders and files will be labelled with this name. "
        "The path to the output folder is set with --out_path.",
    )

    optional.add_argument(
        "-m",
        "--mode",
        dest="mode",
        required=False,
        metavar="MODE",
        help="Specify which BUSCO analysis mode to run.\n"
        "There are three valid modes:\n- geno or genome, for genome assemblies (DNA)\n- tran or "
        "transcriptome, "
        "for transcriptome assemblies (DNA)\n- prot or proteins, for annotated gene sets (protein)",
    )

    optional.add_argument(
        "-l",
        "--lineage_dataset",
        dest="lineage_dataset",
        required=False,
        metavar="LINEAGE",
        help="Specify the name of the BUSCO lineage to be used.",
    )

    optional.add_argument(
        "--augustus",
        dest="use_augustus",
        action="store_true",
        required=False,
        help="Use augustus gene predictor for eukaryote runs",
    )

    optional.add_argument(
        "--augustus_parameters",
        dest="augustus_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        required=False,
        help="Pass additional arguments to Augustus. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--augustus_species",
        dest="augustus_species",
        required=False,
        help="Specify a species for Augustus training.",
    )

    optional.add_argument(
        "--auto-lineage",
        dest="auto-lineage",
        action="store_true",
        required=False,
        help="Run auto-lineage to find optimum lineage path",
    )

    optional.add_argument(
        "--auto-lineage-euk",
        dest="auto-lineage-euk",
        action="store_true",
        required=False,
        help="Run auto-placement just on eukaryote tree to find optimum lineage path",
    )

    optional.add_argument(
        "--auto-lineage-prok",
        dest="auto-lineage-prok",
        action="store_true",
        required=False,
        help="Run auto-lineage just on non-eukaryote trees to find optimum lineage path",
    )

    optional.add_argument(
        "-c",
        "--cpu",
        dest="cpu",
        type=int,
        required=False,
        metavar="N",
        help="Specify the number (N=integer) " "of threads/cores to use.",
    )

    optional.add_argument(
        "--config", dest="config_file", required=False, help="Provide a config file"
    )

    optional.add_argument(
        "--contig_break",
        dest="contig_break",
        type=int,
        required=False,
        metavar="n",
        help="Number of contiguous Ns to signify a break between contigs. Default is n=10.",
    )

    optional.add_argument(
        "--datasets_version",
        dest="datasets_version",
        required=False,
        help="Specify the version of BUSCO datasets, e.g. odb10, odb12 (default odb12)",
    )

    optional.add_argument(
        "--download",
        dest="download",
        required=False,
        type=str,
        metavar="dataset",
        action=DirectDownload,
        help='Download dataset. Possible values are a specific dataset name, "all", "prokaryota", "eukaryota", '
        'or "virus". If used together with other command line arguments, make sure to place this last.',
    )

    optional.add_argument(
        "--download_base_url",
        dest="download_base_url",
        required=False,
        help="Set the url to the remote BUSCO dataset location",
    )

    optional.add_argument(
        "--download_path",
        dest="download_path",
        required=False,
        help="Specify local filepath for storing BUSCO dataset downloads",
    )

    optional.add_argument(
        "-e",
        "--evalue",
        dest="evalue",
        required=False,
        metavar="N",
        type=float,
        help="E-value cutoff for BLAST searches. "
        "Allowed formats, 0.001 or 1e-03 (Default: {:.0e})".format(
            BuscoConfigMain.BLAST_ARGS["evalue"]
        ),
    )

    optional.add_argument(
        "-f",
        "--force",
        action="store_true",
        required=False,
        dest="force",
        help="Force rewriting of existing files. "
        "Must be used when output files with the provided name already exist.",
    )

    optional.add_argument(
        "-h", "--help", action=CleanHelpAction, help="Show this help message and exit"
    )

    optional.add_argument(
        "--limit",
        dest="limit",
        metavar="N",
        required=False,
        type=int,
        help="How many candidate regions (contig or transcript) to consider per BUSCO (default: {})".format(
            str(BuscoConfigMain.BLAST_ARGS["limit"])
        ),
    )

    optional.add_argument(
        "--list-datasets",
        action=ListLineagesAction,
        nargs="?",
        const=None,
        help="Print the list of available BUSCO datasets",
    )

    optional.add_argument(
        "--long",
        action="store_true",
        required=False,
        dest="long",
        help="Optimization Augustus self-training mode (Default: Off); adds considerably to the run "
        "time, but can improve results for some non-model organisms",
    )

    optional.add_argument(
        "--metaeuk",
        dest="use_metaeuk",
        action="store_true",
        required=False,
        help="Use Metaeuk gene predictor",
    )

    optional.add_argument(
        "--metaeuk_parameters",
        dest="metaeuk_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        required=False,
        help="Pass additional arguments to Metaeuk for the first run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--metaeuk_rerun_parameters",
        dest="metaeuk_rerun_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        required=False,
        help="Pass additional arguments to Metaeuk for the second run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    )

    optional.add_argument(
        "--miniprot",
        dest="use_miniprot",
        action="store_true",
        required=False,
        help="Use Miniprot gene predictor",
    )

    optional.add_argument(
        "--skip_bbtools",
        dest="skip_bbtools",
        action="store_true",
        required=False,
        help="Skip BBTools for assembly statistics",
    )

    optional.add_argument(
        "--offline",
        dest="offline",
        action="store_true",
        required=False,
        help="To indicate that BUSCO cannot attempt to download files",
    )

    optional.add_argument(
        "--opt-out-run-stats",
        dest="opt-out-run-stats",
        action="store_true",
        required=False,
        help="Opt out of data collection. Information on the data collected is available in the user guide.",
    )

    optional.add_argument(
        "--out_path",
        dest="out_path",
        required=False,
        metavar="OUTPUT_PATH",
        help="Optional location for results folder, excluding results folder name. "
        "Default is current working directory.",
    )

    optional.add_argument(
        "--plot",
        dest="plot",
        metavar="WORKING_DIRECTORY",
        action=GeneratePlotAction,
        help="Generate a BUSCO summary plot for all short summary files in the given working directory.",
    )

    optional.add_argument(
        "--plot_percentages",
        dest="plot_percentages",
        action="store_true",
        required=False,
        help="Plot the percentages of BUSCOs instead of the number of BUSCOs. To be used as an option with --plot.",
    )

    optional.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        required=False,
        help="Disable the info logs, displays only errors",
        action="store_true",
    )

    optional.add_argument(
        "-r",
        "--restart",
        action="store_true",
        required=False,
        dest="restart",
        help="Continue a run that had already partially completed.",
    )

    optional.add_argument(
        "--scaffold_composition",
        dest="scaffold_composition",
        action="store_true",
        required=False,
        help="Writes ACGTN content per scaffold to a file scaffold_composition.txt",
    )

    optional.add_argument(
        "--tar",
        dest="tar",
        action="store_true",
        required=False,
        help="Compress some subdirectories with many files to save space",
    )

    optional.add_argument(
        "-v",
        "--version",
        action=CleanVersionAction,
        help="Show this version and exit",
        version="BUSCO {}".format(BUSCO.__version__),
    )

    return parser.parse_args(None if len(sys.argv) > 1 else ["--help"])


def main():
    """
    This function runs a BUSCO analysis according to the provided parameters.
    See the help for more details:
    ``BUSCO -h``
    :raises SystemExit: if any errors occur
    """
    params = _parse_args()
    if params.quiet:
        BuscoLogger.quiet = True

    # Execute the GeneratePlotAction if it was triggered
    if hasattr(params, "generate_plot_action"):
        params.generate_plot_action.execute(params)
        sys.exit(0)

    # Execute the DirectDownload action if it was triggered
    if hasattr(params, "direct_download_action"):
        params.direct_download_action.execute(params)
        sys.exit(0)

    busco_run = BuscoMaster(vars(params))
    busco_run.run()


# Entry point
if __name__ == "__main__":
    __spec__ = None
    main()
