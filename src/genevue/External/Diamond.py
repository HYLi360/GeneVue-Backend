from pathlib import Path
from typing import List, Literal

import tomlkit

from genevue.configure import Configure
from genevue.External.CMDBuilder import BatchCMDBuilder

configure = Configure()


class ConfDiamond:
    def __init__(self):
        self.conf_dict = {}

    def load(self):
        pass

    @staticmethod
    def template():
        templ = tomlkit.document()
        templ.add(
            tomlkit.comment(
                "This is a configure template for Diamond. Adjust it as you need!"
            )
        )
        templ.add(
            tomlkit.comment(
                "Note: for params of input and output, please give them as genevue cli arguments, not config!"
            )
        )
        templ.add(tomlkit.nl())

        # diamond makedb: make a reference database in ".dmnd" format.
        makedb = tomlkit.table()
        makedb["taxonmap"] = ""
        makedb["taxonmap"].comment(
            "Gzip file that maps NCBI protein accession numbers to taxon ids.\n"
            "# For Diamond v2.0.7+, Download it at \n"
            "# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz\n"
            "# For others, Download at\n"
            "# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        )
        makedb["taxonnodes"] = ""
        makedb["taxonnodes"].comment(
            "The names.dmp file from the NCBI taxonomy, used to provide taxonomy features.\n"
            "# Download at"
            "# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
        )
        templ.add("makedb", makedb)

        # diamond blastp
        # TODO!

    def save(self):
        pass


def blastp(
    query_seq_path_ls: List[str | Path],
    db_seq_path_ls: List[str | Path],
    res_path: str | Path,
    res_template: str,
    outfmt: str,
    dry_run: bool,
    max_evalue: float | int = 1e-5,
    max_tgt_seqs: int = 0,
    threads: int = 4,
    processes: int = 1,
    verbose: bool = False,
    quiet: bool = True,
    sensitivity: Literal[
        "faster",
        "fast",
        "mid-sensitive",
        "sensitive",
        "more-sensitive",
        "very-sensitive",
        "ultra-sensitive",
    ] = "ultra-sensitive",
):
    cmd = BatchCMDBuilder("diamond", configure.get_program_path("diamond"), "cross")

    query_seq_path_ls = [Path(p) for p in query_seq_path_ls]
    db_seq_path_ls = [Path(p) for p in db_seq_path_ls]
    res_path = Path(res_path)

    cmd.add_flag("blastp")

    cmd.add_substitute_param("--query", query_seq_path_ls)
    cmd.add_substitute_param("--db", db_seq_path_ls)
    cmd.add_substitute_template("--out", "/".join([f"{res_path}", res_template]))
    cmd.add_param("--outfmt", outfmt)
    cmd.add_param("--evalue", max_evalue)

    if max(max_tgt_seqs, 0):
        cmd.add_param("--max_target_seqs", max_tgt_seqs)

    cmd.add_param("--threads", threads)
    cmd.add_flag("--verbose", verbose)
    cmd.add_flag("--quiet", quiet)
    cmd.add_flag(f"--{sensitivity}")

    cmd.run(processes=processes, dry_run=dry_run)
