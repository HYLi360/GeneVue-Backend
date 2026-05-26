from pathlib import Path
from typing import List

from genevue.External.Blaster import BLASTp, MAKEDB
from genevue.External.HMMER3.HMMSearch import HMMSearch
from genevue.Pipelines import Pipeline
from genevue.Sequences.FASTA import FASTA

Pipeline_GeneFamilySearch = Pipeline()


@Pipeline_GeneFamilySearch.node(node_name="HMMsearch", outputs=["hmmsearch_entries"])
def _hmmsearch(
    phmm_probe_path: str, pep_path: str, hmm_res_path: str, hmm_max_evalue: float
):
    hmmsearch_instance = HMMSearch(
        Path(hmm_res_path), Path(phmm_probe_path), Path(pep_path)
    )
    hmmsearch_instance.run()
    return hmmsearch_instance.result_entries_filter(hmm_max_evalue)


@Pipeline_GeneFamilySearch.node(node_name="MakeBLASTDB")
def _makeblastdb(pep_path: str, blastp_pdb_path: str):
    makedb_instance = MAKEDB("diamond", Path(pep_path), Path(blastp_pdb_path))
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
    cross_set = set(hmmsearch_entries) & set(blastp_entries)
    with open("cross.txt", "w") as f:
        for cross in cross_set:
            f.write(f"{cross}\n")

    pep_fasta = FASTA(Path(pep_path))
    pep_fasta.filter(Path("cross.txt"), "plain", Path(record_output_path))
