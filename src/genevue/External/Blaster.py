from typing import Literal
from pathlib import Path
import subprocess
import shlex

from multiprocessing import cpu_count
from dataclasses import dataclass

from genevue.configure import Configure
from genevue import setup_rich_logger, console
from genevue.utils.parse import blast6reader
from genevue.External.CMDBuilder import CMDBuilder

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
        db_path: Path,
        verbose: bool = False,
        threads: int | Literal["max", "single"] = 4,
    ):
        self.program_name = method
        self.program_path = configure.get_program_path(method)
        self.dbseqs_path = dbseqs_path
        self.db_path = db_path
        self.verbose = verbose

        if isinstance(threads, int):
            self.threads = threads
        else:
            self.threads = {"max": CPU_COUNT, "single": 1}.get(threads)

        self.cmdbuilder = CMDBuilder(self.program_name, self.program_path)

        self.buildcmd()

    def buildcmd(self):
        match self.program_name:
            case "makeblastdb":
                pass
            case "diamond":
                self.cmdbuilder.add_flag("makedb")
                self.cmdbuilder.add_param("--in", self.dbseqs_path)
                self.cmdbuilder.add_param("--db", self.db_path)
                self.cmdbuilder.add_param("--threads", self.threads)
                self.cmdbuilder.add_flag("--quiet")
                self.cmdbuilder.add_flag("--verbose", self.verbose)
                return
        # this branch should not be reached
        logger.warning(f"Not identified MAKEDB method: {self.program_name}")

    def run(self, dry_run: bool = False):
        self.cmdbuilder.run(dry_run=dry_run)


class BLASTp:
    def __init__(
        self,
        method: Literal["NCBI BLAST+", "DIAMOND", "MMseqs2"],
        query_seqs_path: str | Path,
        db_path: str | Path,
        res_path: str | Path,
    ) -> None:
        self.method = method
        self.query_seqs_path: Path = (
            Path(query_seqs_path).resolve()
            if isinstance(query_seqs_path, str)
            else query_seqs_path
        )
        self.db_path: Path = (
            Path(db_path).resolve() if isinstance(db_path, str) else db_path
        )
        self.res_path: Path = (
            Path(res_path).resolve() if isinstance(res_path, str) else res_path
        )
        self.res_format: Literal["pairwise", "xml", "tsv"] | int = "tsv"
        self.max_evalue: float = 1e-5
        self.threads: Literal["auto", "single"] | int = "auto"
        self.repeat_masking: bool = True
        self.min_bitscore: float = 0
        self.max_tgt_seqs: int = 10000
        self.min_identify_percent: float = 0

        self.diamond_sensitivity: Literal[
            "faster",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ] = "ultra-sensitive"

        # Internal label for genevue identifying
        self._input_label = ["seqs", "blastdb"]
        self._output_label = "blastres"

        # for res_format (outfmt)
        # if using literal, turn it into int
        if not isinstance(self.res_format, int):
            self.res_format = {"pairwise": 0, "xml": 5, "tsv": 6}.get(
                self.res_format, self.res_format
            )

        # for threads
        # if using "auto", transform it to maximum thread number
        # using "single", transform it to single thread
        if not isinstance(self.threads, int):
            self.threads = {"auto": CPU_COUNT, "single": 1}.get(
                self.threads, self.threads
            )

    @property
    def command(self):
        match self.method:
            case "NCBI BLAST+":
                return []
            case "DIAMOND":
                return [
                    "diamond",
                    "blastp",
                    "--query",
                    f"{self.query_seqs_path}",
                    "--db",
                    f"{self.db_path}",
                    "--out",
                    f"{self.res_path}",
                    "--evalue",
                    f"{self.max_evalue}",
                    "--threads",
                    f"{int(self.threads)}",
                    "--max-target-seqs",
                    f"{self.max_tgt_seqs}",
                    "--min-score",
                    f"{self.min_bitscore}",
                    "--outfmt",
                    f"{self.res_format}",
                    f"--{self.diamond_sensitivity}",
                ]
            case "MMseqs2":
                return []
        # this branch should not be reached
        print(f"Not identified blastp method: {self.method}")
        return []

    def run(self):
        print(f"{self.method} BLASTp started.")
        print(f"cmd: {shlex.join(self.command)}")
        subprocess.run(self.command)

    def parse_blast6(self, hi_evalue: float, lo_bitscore: float):
        return blast6reader(self.res_path, hi_evalue, lo_bitscore)

    def parse_blast5(self):
        pass
