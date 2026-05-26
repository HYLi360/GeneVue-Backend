""" """

""" """

import re
from pathlib import Path
from typing import Optional

import pandas as pd

from genevue import setup_rich_logger, console
from genevue.External.CMDBuilder import CMDBuilder
from genevue.Pipelines import ProcessNode
from genevue.configure import Configure

logger = setup_rich_logger(__name__, console)
configure = Configure()


class HMMSearch(ProcessNode):
    def __init__(
        self,
        res_path: Path,
        phmm_probe_path: Optional[Path] = None,
        seqs_path: Optional[Path] = None,
    ) -> None:
        self.phmm_probe_path = phmm_probe_path
        self.seqs_path = seqs_path
        self.res_path = res_path
        self.command = CMDBuilder("hmmsearch", configure.get_program_path("hmmsearch"))
        self.result_info = None
        self.result_table = None
        self.result_entries = None

        self.buildcmd()
        self.run()

    def buildcmd(self):
        self.command.add_param("-o", self.res_path)

    def run(self):
        self.command.run(self.phmm_probe_path, self.seqs_path)
        self.parse_result()

    def parse_result(self):
        regexls = [
            r"^# HMMER ([\d\.]+).+$",
            r"^# target sequence database:\s+(.+)$",
            r"^Query:\s+(.+)$",
            r"^([\d\.\-e]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.\-e]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+\d+\s+[^\s]+\s+.+)$",
            r"^Target sequences:\s+(\d+).+$",
            r"^Domain search space  \(domZ\):\s+(\d+).+$",
        ]

        res_containor = []
        with open(self.res_path, "r") as f:
            for line in f:
                line = line.strip()
                for regex in regexls:
                    group = re.search(regex, line)
                    if group:
                        res_containor.append(group.group(1))

            self.result_info = dict(
                zip(
                    [
                        "version",
                        "target",
                        "query",
                        "target sequences",
                        "passed sequences",
                    ],
                    [
                        res_containor[0],
                        res_containor[1],
                        res_containor[2],
                        res_containor[-2],
                        res_containor[-1],
                    ],
                )
            )

            res_containor = res_containor[3:-2]
            pre_table = []
            for line in res_containor:
                items = line.split(maxsplit=9)
                items[0:7] = [float(value) for value in items[0:7]]
                pre_table.append(items)

            self.result_table = pd.DataFrame(
                data=pre_table,
                columns=[
                    ("full", "E-value"),
                    ("full", "score_full"),
                    ("full", "bias"),
                    ("best", "E-value"),
                    ("best", "score"),
                    ("best", "bias"),
                    ("domains", "exp"),
                    ("domains", "N"),
                    "Sequence",
                    "Description",
                ],
            )

            self.result_entries = self.result_table["Sequence"].tolist()

    def result_entries_filter(self, max_evalue: Optional[float] = None):
        if max_evalue:
            return self.result_table[
                self.result_table[("full", "E-value")] <= max_evalue
            ]["Sequence"].tolist()
        return self.result_entries
