from pathlib import Path
from typing import List

import Bio.SeqIO

from genevue.external.Blaster import BLASTp, MAKEDB
from genevue.external.HMMER3 import HMMSEARCH
from genevue.Pipelines import Pipeline
from genevue.utils.parse import records_filter

Pipeline_GeneFamilySearch = Pipeline()


@Pipeline_GeneFamilySearch.node(node_name="HMMsearch", outputs=["hmmsearch_entries"])
def _hmmsearch(
    phmm_probe_path: str, pep_path: str, hmm_res_path: str, hmm_max_evalue: float
):
    hmmsearch_instance = HMMSEARCH(phmm_probe_path, pep_path, hmm_res_path)
    hmmsearch_instance.run()
    return hmmsearch_instance.result_entries_filter(hmm_max_evalue)


@Pipeline_GeneFamilySearch.node(node_name="MakeBLASTDB")
def _makeblastdb(pep_path: str, blastp_pdb_path: str):
    makedb_instance = MAKEDB("DIAMOND", pep_path, blastp_pdb_path)
    makedb_instance.run()


@Pipeline_GeneFamilySearch.node(
    node_name="BLASTp", outputs=["blastp_entries"], depends=["MakeBLASTDB"]
)
def _blastp(
    blastp_probe_path: str,
    blastp_pdb_path: str,
    blastp_res_path: str,
    blastp_max_evalue: float,
    blastp_min_bitscore: float,
):
    blastp_instance = BLASTp(
        "DIAMOND", blastp_probe_path, blastp_pdb_path, blastp_res_path
    )
    blastp_instance.run()
    return blastp_instance.parse_blast6(blastp_max_evalue, blastp_min_bitscore)[
        "gene2"
    ].tolist()


@Pipeline_GeneFamilySearch.node(
    node_name="extract", outputs=["records"], depends=["HMMsearch", "BLASTp"]
)
def _extract(
    pep_path: str,
    hmmsearch_entries: List[str],
    blastp_entries: List[str],
    record_output_path: str,
):
    cross = set(hmmsearch_entries) & set(blastp_entries)
    with open(Path(record_output_path).resolve(), "w+") as f:
        Bio.SeqIO.write(
            records_filter(Path(pep_path).resolve(), list(cross)).values(), f, "fasta"
        )
    return list(records_filter(Path(pep_path).resolve(), list(cross)).values())
