#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

"""
UniPort download interface based on Proteins API.

Andrew Nightingale, Ricardo Antunes, Emanuele Alpi, Borisas Bursteinas, Leonardo Gonzales, Wudong Liu, Jie Luo,
Guoying Qi, Edd Turner, Maria Martin, The Proteins API: accessing key integrated protein and genome information,
Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W539-W544. https://doi.org/10.1093/nar/gkx237
"""
from types import NoneType
from typing import Literal

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import requests

from genevue import console


class Protein:
    """
    Protein Service by Protein API.

    key notes on some parameters:

    - `offset` and `size`(int). When search results are numerous, you can set the `size` parameter to control the
    number of results returned. The `offset` parameter specifies the starting position within the search results (0-based).

    - `accession`(List), UniProt accessions. it is (they are) the unique identifier(s) for protein entries in
    UniPortKB, consisting of 6 or 10 characters (e.g., A2BC19, P12345, A0A023GPI8). Up to 100.

    - `reviewed`(Bool), reviewed label of entries.

    - `isoform`(String), 0 for no isoform, 1 for isoform only, and 2 for both canonical and isoform. for semenic
    option, use "no", "only" or "both".

    - `goterms`(String), Gene Ontology term (name only. GO accession is unacceptable).

    - `keywords`(String), word in UniPort Keywords. UniPort Keywords is a controlled vocabulary table,
    like Gene Ontology or NLM MeSH. https://www.uniprot.org/keywords/

    - `ec`(List), EC number(s). EC (Enzyme Commission) number is a standard for coding every enzymes. Up to 20.

    - `gene`/`exact_gene`/`protein`(List), Gene/exact-gene/protein name(s) in UniProt. Up to 20.

    - `organism`(String)/`taxid`(List), Organism name/taxon ID. `organism` only receives one organism name,
    but `taxid` can receives up to 20.

    - `pubmed`(String), UniProt reference PubMed ID.

    - `seqLength`(String), value (123) or range (123-456) for query sequence(s).

    - `md5`(String), sequence md5 value.

    The above parameters can be entered simultaneously. However, please note that these parameters are linked by an
    “AND” relationship: search results must satisfy all parameters.
    """

    def __init__(self):
        self._baseurl = "https://www.ebi.ac.uk/proteins/api/proteins"
        self._headers = {"Accept": "text/x-fasta"}
        self.request_url = None
        self.result = None
        self.seqrecords = {}
        console.log("remote.ProteinAPI.Protein called.")

    def exactly_search(
        self,
        offset: int = 0,
        size: int = 100,
        accession: list | None = None,
        reviewed: bool | None = None,
        isoform: Literal["0", "1", "2", "no", "only", "both"] | None = None,
        goterms: str | None = None,
        keywords: str | None = None,
        ec: list | None = None,
        gene: list | None = None,
        exact_gene: list | None = None,
        protein: list | None = None,
        organism: str | None = None,
        taxid: list | None = None,
        pubmed: list | None = None,
        seqLength: str | None = None,
        md5: str | None = None,
    ) -> None:
        # re-format params
        paramls = [
            f"offset={offset}",
            f"size={size}",
            (
                f"accession={"%2C".join(accession)}"
                if not isinstance(accession, NoneType)
                else None
            ),
            f"reviewed={reviewed}" if not isinstance(reviewed, NoneType) else None,
            (
                f"isoform={({"no": "0", "only": "1", "both": "2"}.get(isoform, isoform))}"
                if not isinstance(isoform, NoneType)
                else None
            ),
            (
                f"goterms={goterms}".replace(" ", "%20")
                if not isinstance(goterms, NoneType)
                else None
            ),
            (
                f"keywords={keywords}".replace(" ", "%20")
                if not isinstance(keywords, NoneType)
                else None
            ),
            f"ec={"%2C".join(ec)}" if not isinstance(ec, NoneType) else None,
            f"gene={"%2C".join(gene)}" if not isinstance(gene, NoneType) else None,
            (
                f"exact_gene={"%2C".join(exact_gene)}"
                if not isinstance(exact_gene, NoneType)
                else None
            ),
            (
                f"protein={"%2C".join(protein)}"
                if not isinstance(protein, NoneType)
                else None
            ),
            (
                f"organism={"%2C".join(organism)}".replace(" ", "%20")
                if not isinstance(organism, NoneType)
                else None
            ),
            f"taxid={"%2C".join(taxid)}" if not isinstance(taxid, NoneType) else None,
            (
                f"pubmed={"%2C".join(pubmed)}"
                if not isinstance(pubmed, NoneType)
                else None
            ),
            f"seqLength={seqLength}" if not isinstance(seqLength, NoneType) else None,
            f"md5={md5}" if not isinstance(md5, NoneType) else None,
        ]
        # remove None params
        paramls: list[str] = [param for param in paramls if param is not None]
        # request url
        self.request_url = self._baseurl + "?" + "&".join(paramls)
        with console.status("Downloading data from EMBL-EBI Proteins API service."):
            console.log(f"Request URL: {self.request_url}")
            console.log(f"Cite the Proteins API!")
            console.log(
                "The proteins API: accessing key integrated protein and genome information (doi:10.1093/nar/gkx237)"
            )

            request = requests.get(self.request_url, headers=self._headers)
            if not request.ok:
                print(request.json())
                request.raise_for_status()
        self.result = request.json()
        self._transform()

    def batch_fetch(self, accessions: list):
        offset = 0
        while len(accessions) - offset >= 100:
            self.exactly_search(accession=accessions[offset : offset + 100])
            offset += 100
        if len(accessions) - offset != 0:
            self.exactly_search(accession=accessions[offset:])

    def _transform(self):
        """Transform the search results to process-able SeqRecord dictionary."""
        pass

    def _parse_query(self, file_path: str):
        pass
