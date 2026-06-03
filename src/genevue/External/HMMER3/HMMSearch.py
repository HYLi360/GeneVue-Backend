""" """

from collections import defaultdict


import re
from pathlib import Path
from typing import Optional

import pandas as pd

from genevue import setup_rich_logger, console
from genevue.External.CMDBuilder import CMDBuilder, BatchCMDBuilder
from genevue.Pipelines import ProcessNode
from genevue.configure import Configure
from genevue.Tools.filesystem import iter_path

from Bio import SeqIO

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


class HMMSearchResultDomtblout:
    def __init__(self, result_path: Path, seq_output_path: Path):
        self.path = result_path
        self.seq_output_path = seq_output_path
        self._parse()

    def _parse(self):
        res = defaultdict(list)
        with open(self.path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    # for these 23 columns:
                    #      0 = target sequence name (may duplicated)
                    #      1 = target accession (may duplicated)
                    #      2 = target length
                    #      3 = query hmm name
                    #      4 = query accession
                    #      5 = query length
                    #      6 = full-sequence evalue
                    #      7 =               score
                    #      8 =               bias
                    #      9 = domain order
                    #     10 = domains count
                    #     11 = domain evalue
                    #     12 =        i-evalue
                    #     13 =        score
                    #     14 =        bias
                    # 15, 16 = align range on hmm coord
                    # 17, 18 = align range on query coord
                    # 19, 20 = align range on query coord (with envolope)
                    #     21 = accuraty
                    #     22 = query sequence description
                    resls = line.strip().split(maxsplit=22)
                    res[resls[0]].append(
                        {
                            "hmm name": resls[3],
                            "hmm accession": resls[4],
                            "hmm range": (int(resls[15]), int(resls[16])),
                            "target range": (int(resls[19]), int(resls[20])),
                            "domain order": int(resls[9]),
                            "domain evalue": float(resls[11]),
                            "domain ievalue": float(resls[12]),
                            "domain score": float(resls[13]),
                            "domain bias": float(resls[14]),
                            "accuraty": float(resls[21]),
                            "description": f"hmmsearch_{resls[3]}@{resls[15]}-{resls[16]}_{resls[0]}@{resls[19]}-{resls[20]}"
                            f"_E{resls[11]}_S{resls[13]}_A{resls[21]}",
                        }
                    )
        self.result = dict(res)


class BatchHMMSearch:
    def __init__(
        self,
        seq_dir: Path,
        res_dir: Path,
        extract_dir: Path,
        phmm_probe_path: Path,
    ) -> None:
        self.phmm_probe_path = phmm_probe_path
        self.seq_dir = seq_dir
        self.res_dir = res_dir
        self.extract_dir = extract_dir
        self.command = BatchCMDBuilder(
            "hmmsearch",
            configure.get_program_path("hmmsearch"),
            substitute_method="paired",
        )
        self.seq_file_ls = iter_path(target_path=self.seq_dir)
        self.out_path_ls = [
            res_dir / f"{p.stem}_hmmsearch_res.txt" for p in self.seq_file_ls
        ]
        self.command.add_substitute_param("--domtblout", self.out_path_ls)
        self.command.add_flag(f"{phmm_probe_path}")
        self.command.add_substitute_param(
            param_name=None, substitute_list=self.seq_file_ls
        )

    def run(self, dry_run=True):
        self.command.run(dry_run=dry_run)
        self.extract()

    def extract(self):
        logger.info("Start sequence extracting.")
        for seq, res in zip(self.seq_file_ls, self.out_path_ls):
            logger.info(f"Extracting from {seq}")
            resdict = defaultdict(list)
            with open(res) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    else:
                        # for these 23 columns:
                        #      0 = target sequence name (may duplicated)
                        #      1 = target accession (may duplicated)
                        #      2 = target length
                        #      3 = query hmm name
                        #      4 = query accession
                        #      5 = query length
                        #      6 = full-sequence evalue
                        #      7 =               score
                        #      8 =               bias
                        #      9 = domain order
                        #     10 = domains count
                        #     11 = domain evalue
                        #     12 =        i-evalue
                        #     13 =        score
                        #     14 =        bias
                        # 15, 16 = align range on hmm coord
                        # 17, 18 = align range on query coord
                        # 19, 20 = align range on query coord (with envolope)
                        #     21 = accuraty
                        #     22 = query sequence description
                        resls = line.strip().split(maxsplit=22)
                        resdict[resls[0]].append(
                            f"hmmsearch_{resls[3]}@{resls[15]}-{resls[16]}_{resls[0]}@{resls[19]}-{resls[20]}"
                            f"_E{resls[11]}_S{resls[13]}_A{resls[21]}",
                        )
            resdict = dict(resdict)
            with open(f"{self.extract_dir}/{seq.stem}.faa", "w") as f:
                for record in SeqIO.parse(seq, "fasta"):
                    if record.id in resdict:
                        record.description = (
                            f"{record.description} {" ".join(resdict[record.id])}"
                        )
                        SeqIO.write(record, f, "fasta")
