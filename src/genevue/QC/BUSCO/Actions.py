# coding: utf-8
"""
Actions.py

Custom command line actions.

Author(s): Matthew Berkeley

Copyright (c) 2015-2025, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

import argparse
import os
import sys

from genevue.QC.BUSCO.BuscoConfig import PseudoConfigForDownload, PseudoConfigForPlot
from genevue.QC.BUSCO.BuscoDownloadManager import BuscoDownloadManager
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
from genevue.QC.BUSCO.ConfigManager import BuscoConfigManager


class ListLineagesAction(argparse.Action):

    logger = BuscoLogger.get_logger(__name__)

    def __init__(
        self, option_strings, dest, nargs="?", default="==SUPPRESS==", **kwargs
    ):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            self.config_manager = BuscoConfigManager(namespace.__dict__)
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)

        self.config = PseudoConfigForDownload(
            self.config_manager.config_file, self.config_manager.params
        )
        try:
            self.config.load()
            self.print_lineages(values)
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)
        finally:
            try:
                os.remove("busco_{}.log".format(BuscoLogger.pid))
            except FileNotFoundError:
                pass
            parser.exit()

    def print_lineages(self, values):
        if values is not None:
            if values == "odb10":
                lineages_list_version = "odb10"
            elif values == "odb12":
                lineages_list_version = "odb12"
            else:
                type(self).logger.error(
                    "{} is not a recognized option".format(values[0])
                )
                return
        else:
            lineages_list_version = "odb12"

        if self.config.update:
            lineages_list_file = self.download_lineages_list(lineages_list_version)
        else:
            lineages_list_file = self.config.existing_downloads[0]
        with open(lineages_list_file, "r") as f:
            print("".join(f.readlines()))

    def download_lineages_list(self, version):
        if version == "odb10":
            lineages_list_file = self.config.downloader.get(
                "lineages_list_odb10.txt", "information"
            )
        elif version == "odb12":
            lineages_list_file = self.config.downloader.get(
                "lineages_list.txt", "information"
            )
        else:
            type(self).logger.error("{} is not a recognized option".format(version))
            return
        return lineages_list_file


class DirectDownload(argparse.Action):
    logger = BuscoLogger.get_logger(__name__)

    def __init__(
        self, option_strings, dest, nargs="*", default="==SUPPRESS==", **kwargs
    ):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        # Store the action and parameters in the namespace for later execution
        setattr(namespace, "direct_download_action", self)
        setattr(namespace, "download_values", values)

    def execute(self, namespace):
        try:
            self.config_manager = BuscoConfigManager(vars(namespace))
            self.config = PseudoConfigForDownload(
                self.config_manager.config_file, self.config_manager.params
            )
            self.config.load()
            bdm = BuscoDownloadManager(self.config)
            self.download_datasets(bdm, namespace.download_values)
        except SystemExit as se:
            type(self).logger.error(se)
            raise SystemExit(1)
        finally:
            try:
                os.remove("busco_{}.log".format(BuscoLogger.pid))
            except FileNotFoundError:
                pass

    def download_datasets(self, bdm, values):
        for item in values:
            self.get(bdm, item)

    def get(self, bdm, item):
        if item == "all" or item == "all_odb12":
            files_to_get = [item for item in bdm.version_files.index if "odb12" in item]
        elif item == "prokaryota" or item == "all_prokaryota_odb12":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if (bdm.version_files.loc[item]["domain"] == "Prokaryota")
                and ("odb12" in item)
            ]
        elif item == "eukaryota" or item == "all_eukaryota_odb12":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if (bdm.version_files.loc[item]["domain"] == "Eukaryota")
                and ("odb12" in item)
            ]
        elif item == "virus":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if bdm.version_files.loc[item]["domain"] == "Virus"
            ]
        elif item == "all_odb10":
            files_to_get = [item for item in bdm.version_files.index if "odb10" in item]
        elif item == "all_prokaryota_odb10":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if (bdm.version_files.loc[item]["domain"] == "Prokaryota")
                and ("odb10" in item)
            ]
        elif item == "all_eukaryota_odb10":
            files_to_get = [
                item
                for item in bdm.version_files.index
                if (bdm.version_files.loc[item]["domain"] == "Eukaryota")
                and ("odb10" in item)
            ]
        else:
            try:
                if isinstance(item, str) and item in bdm.version_files.index:
                    files_to_get = [item]
                else:
                    raise KeyError
            except KeyError:
                type(self).logger.error("{} is not a recognized option".format(item))
                files_to_get = []

        filetypes = [bdm.version_files.loc[f]["type"] for f in files_to_get]
        for f, filename in enumerate(files_to_get):
            bdm.get(filename, filetypes[f])


class CleanHelpAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=0, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        try:
            os.remove("busco_{}.log".format(BuscoLogger.pid))
        except OSError:
            pass
        parser.exit()


class CleanVersionAction(argparse.Action):
    def __init__(
        self,
        option_strings,
        version=None,
        dest="==SUPPRESS==",
        nargs=0,
        default="==SUPPRESS==",
        **kwargs,
    ):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)
        self.version = version

    def __call__(self, parser, namespace, values, option_string=None):
        version = self.version
        if version is None:
            version = "unknown"
        formatter = parser._get_formatter()
        formatter.add_text(version)
        parser._print_message(formatter.format_help(), sys.stdout)
        try:
            os.remove("busco_{}.log".format(BuscoLogger.pid))
        except OSError:
            pass
        parser.exit()


class GeneratePlotAction(argparse.Action):

    def __init__(self, option_strings, dest, nargs=1, default="==SUPPRESS==", **kwargs):
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        # Store the parameters in the namespace for later execution
        setattr(namespace, "generate_plot_action", self)
        setattr(namespace, "plot_values", values)

    def execute(self, namespace):
        try:
            from genevue.QC.BUSCO.BuscoPlot import Plot

            self.config_manager = BuscoConfigManager(vars(namespace))
            self.config = PseudoConfigForPlot(
                self.config_manager.config_file, self.config_manager.params
            )
            self.config.load()
            plot = Plot(self.config, namespace.plot_values[0])
            plot.load_data()
            plot.generate_plot()
        except SystemExit as se:
            BuscoLogger.get_logger(__name__).error(se)
            raise SystemExit(1)
        finally:
            try:
                os.remove("busco_{}.log".format(BuscoLogger.pid))
            except FileNotFoundError:
                pass
