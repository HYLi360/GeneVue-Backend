import os
from typing import Literal, Optional
from pathlib import Path

from multiprocessing import cpu_count
from dataclasses import dataclass

from genevue.Tools import iter_path
from genevue.configure import Configure
from genevue import setup_rich_logger, console
from genevue.External.CMDBuilder import BatchCMDBuilder

CPU_COUNT = cpu_count()

logger = setup_rich_logger(__name__, console)

configure = Configure()


@dataclass(
    kw_only=True,
)
class MAKEDB:
    def __init__(
        self,
        method: Literal["makeblastdb", "diamond"],
        dbseqs_path: Path,
        out_db_path: Path,
        db_type: Literal["nucl", "prot"],
        verbose: bool = False,
        threads: int | Literal["max", "single"] = 4,
        processes: int = 1,
    ) -> None:
        self.program_name = method
        self.program_path = configure.get_program_path(method)
        self.dbseqs_path = dbseqs_path
        self.out_db_path = out_db_path
        self.db_type = db_type
        self.verbose = verbose

        if isinstance(threads, int):
            self.threads = threads
        else:
            self.threads = {"max": CPU_COUNT, "single": 1}.get(threads)

        self.processes = processes

        # Internal label for genevue identifying
        self._input_label = ["seqs"]
        self._output_label = f"blastdb_{db_type}"

        self._ready = False

        self.cmd = BatchCMDBuilder(self.program_name, self.program_path)

        self._buildcmd()

    def _buildcmd(self) -> None:
        match self.program_name:
            case "makeblastdb":
                return
            case "diamond":
                if self.db_type != "prot":
                    logger.error("Program diamond can only make protein database")
                    return
                self.cmd.add_flag("makedb")

                if self.dbseqs_path.is_dir():
                    in_file_list = iter_path(self.dbseqs_path)
                else:
                    in_file_list = [self.dbseqs_path]

                self.cmd.add_substitute_param("--in", in_file_list)
                self.cmd.add_substitute_template(
                    "--db", "/".join([f"{self.out_db_path}", "{0.stem}.dmnd"])
                )
                self.cmd.add_param("--threads", self.threads)
                self.cmd.add_flag("--quiet")
                self.cmd.add_flag("--verbose", self.verbose)
                self._ready = True
                return

    def run(self, dry_run: bool = False) -> None:
        if self._ready:
            self.cmd.run(dry_run=dry_run, processes=self.processes)
        else:
            logger.error("Prevent execution due to expection happened before")

    def dry_run(self):
        """
        Simplified of self.run(dry_run=True).
        """
        self.run(dry_run=True)


class BLASTp:
    def __init__(
        self,
        method: Literal["blastp", "diamond"],
        query_seqs_path: Path,
        db_path: Path,
        res_path: Path,
        max_evalue: float | int = 1e-5,
        min_bitscore: float | int = 0,
        threads: int | Literal["max", "single"] = 4,
        repeat_masking: bool = True,
        max_tgt_seqs: int = 0,
        output_format: Literal["xml", "tsv"] | int = "tsv",
        verbose: bool = False,
        quiet: bool = True,
        processes: int = 4,
        other_params: Optional[list] = None,
        diamond_sensitivity: Literal[
            "faster",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ] = "ultra-sensitive",
    ) -> None:
        self.program_name = method
        self.program_path = configure.get_program_path(method)
        self.method = method
        self.query_seqs_path = query_seqs_path
        self.db_path = db_path
        self.res_path = res_path
        self.max_evalue: float | int = max_evalue
        self.min_bitscore: float | int = min_bitscore

        if isinstance(threads, int):
            self.threads = threads
        else:
            self.threads = {"max": CPU_COUNT, "single": 1}.get(threads)

        self.repeat_masking = repeat_masking
        self.max_tgt_seqs = max_tgt_seqs

        if isinstance(output_format, int):
            self.output_format = output_format
        else:
            self.output_format = {"xml": 5, "tsv": 6}.get(output_format)

        self.verbose = verbose
        self.quiet = quiet
        self.processes = processes
        self.other_params = other_params

        self.diamond_sensitivity = diamond_sensitivity

        # Internal label for genevue identifying
        self._input_label = ["seqs", "blastdb"]
        self._output_label = "blastres"

        self.cmd = BatchCMDBuilder(
            self.program_name, self.program_path, substitute_method="cross"
        )

        self._buildcmd()

    def _buildcmd(self):
        match self.method:
            case "blastp":
                self.cmd.add_param("-query", self.query_seqs_path)
                self.cmd.add_param("-db", self.db_path)
                self.cmd.add_param("-out", self.res_path)
                self.cmd.add_param("-outfmt", self.output_format)
                self.cmd.add_param("-evalue", self.max_evalue)
                self.cmd.add_param("-max_target_seqs", self.max_tgt_seqs)
                self.cmd.add_param("-threads", self.threads)

                # custom params
                if self.other_params:
                    self.cmd.add_custom_params(self.other_params)
                return
            case "diamond":
                self.cmd.add_flag("blastp")

                in_seq_file_list = iter_path(self.query_seqs_path)
                in_db_file_list = iter_path(self.db_path)

                self.cmd.add_substitute_param("--query", in_seq_file_list)
                self.cmd.add_substitute_param("--db", in_db_file_list)
                self.cmd.add_substitute_template(
                    "--out", "/".join([f"{self.res_path}", "{0.stem}_{1.stem}.out"])
                )
                self.cmd.add_param("--outfmt", self.output_format)
                self.cmd.add_param("--evalue", self.max_evalue)

                if self.max_tgt_seqs > 0:
                    self.cmd.add_param("--max_target_seqs", self.max_tgt_seqs)

                self.cmd.add_param("--threads", self.threads)
                self.cmd.add_flag("--verbose", self.verbose)
                self.cmd.add_flag("--quiet", self.quiet)
                self.cmd.add_flag(f"--{self.diamond_sensitivity}")

                # custom params
                if self.other_params:
                    self.cmd.add_custom_params(self.other_params)
                return

    def run(self, dry_run: bool = False):
        logger.info(f"{self.method} BLASTp started.")
        self.cmd.run(dry_run=dry_run, processes=self.processes)

    def dry_run(self):
        """
        Simplified of self.run(dry_run=True).
        """
        self.run(dry_run=True)
