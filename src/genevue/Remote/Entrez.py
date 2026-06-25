#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

"""
A simple tool to execute single (or batch) search on NCBI Entrez(R). NCBI Entrez(R) is a service provided by NCBI, which
make users easily to search cross-database.

When using this tool, please ensure you have a stable internet connection, provide your email address (by `genevue config
email`) and API key (`genevue config api-key NCBI [your-api-key-here]`, optional. NCBI account needed), and agree to
Entrez's Usage Guidelines (see at <https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_>).
"""

from pathlib import Path
from typing import Optional, List, Dict

import pandas as pd
import typer
from Bio import Entrez as BioEntrez
from rich.text import Text

from genevue import (
    setup_rich_logger,
    console,
    AllFieldsEmptyError,
)
from genevue.Utils.Parse import species_italic_name
from genevue.configure import Configure

logger = setup_rich_logger(__name__, console)
configure = Configure()
BioEntrez.email = configure.email
app_entrez = typer.Typer()


class GVEUtilsTaxonomy:
    def __init__(
        self,
        species_id: int | str = "",
        species_name: str = "",
    ):
        self.species_id: str = str(species_id)
        self.species_name: str = species_name

        # Slot for a search result
        self.genetic_code: Optional[int] = None
        self.genetic_code_mt: Optional[int] = None

        self.lineage_dict: Dict[str, str] = {}

        self._isok = False

    def get(self):
        subterm_list = []
        for subterm_content, subterm_flag in [
            (self.species_id, "[UID]"),
            (self.species_name, "[SCIN]"),
        ]:
            if subterm_content:
                subterm_list.append("".join([subterm_content, subterm_flag]))

        if not subterm_list:
            logger.exception(
                AllFieldsEmptyError("species_id", "species_name"),
            )

        # Get ID
        stream = BioEntrez.esearch(db="taxonomy", term=" AND ".join(subterm_list))
        record = BioEntrez.read(stream)
        if not record["IdList"]:
            # try to get the true spell
            stream = BioEntrez.espell(term=self.species_name)
            correct_spell = BioEntrez.read(stream)["CorrectedQuery"]
            if not self.species_name or not correct_spell:
                logger.error(
                    "I tried searching but found nothing. Maybe you should check your input?"
                )
                return
            else:
                logger.error(
                    "I tried searching but found nothing. Maybe you spelled wrong, "
                    f"mixed [red]{self.species_name.lower()}[/red] with [green]{correct_spell}[/green] ?"
                )
                return

        stream = BioEntrez.efetch(db="taxonomy", id=record["IdList"][0])
        record = BioEntrez.read(stream)[0]
        self.species_id = record["TaxId"]
        self.species_name = record["ScientificName"]
        self.genetic_code = record["GeneticCode"]["GCId"]
        self.genetic_code_mt = record["MitoGeneticCode"]["MGCId"]

        if record["Rank"] != "species":
            logger.error(f"Result '{record["ScientificName"]}' is not a species")
            return

        for rank in record["LineageEx"]:
            match rank["Rank"]:
                case "genus":
                    self.lineage_dict["genus"] = rank["ScientificName"]
                case "family":
                    self.lineage_dict["family"] = rank["ScientificName"]
                case "order":
                    self.lineage_dict["order"] = rank["ScientificName"]
                case "class":
                    self.lineage_dict["class"] = rank["ScientificName"]

        self._isok = True

    def show(self):
        if not self.lineage_dict:
            self.get()

        if self._isok:
            console.print(
                species_italic_name(self.species_name), f"[{self.species_id}]"
            )
            console.print(f"Class:  {self.lineage_dict['class']}")
            console.print(f"Order:  {self.lineage_dict['order']}")
            console.print(f"Family: {self.lineage_dict['family']}")
            console.print(
                Text("Genus: "), Text(self.lineage_dict["genus"], style="italic")
            )
            console.print(f"Nucl-Genetic Code Type: {self.genetic_code}")
            console.print(f"MT - Genetic Code Type: {self.genetic_code_mt}")

    @property
    def lineage(self) -> Dict[str, str]:
        if not self.lineage_dict:
            self.get()
        if self._isok:
            return {
                "Class": self.lineage_dict["class"],
                "Order": self.lineage_dict["order"],
                "Family": self.lineage_dict["family"],
                "Genus": self.lineage_dict["genus"],
                "Species": self.lineage_dict["Species"],
            }
        return {}


class GVEUtilsBatchTaxonomy:
    def __init__(
        self,
        species_name_ls: Optional[List[str]] = None,
        species_id_ls: Optional[List[str]] = None,
    ):
        self.species_name_ls = species_name_ls if species_name_ls is not None else []
        self.species_id_ls = species_id_ls if species_id_ls is not None else []

        # Slot for a search result
        self.lineage_list: List[Dict[str, str]] = []
        self._isok = False

    def get(self):
        logger.info(
            f"Get {len(self.species_name_ls) + len(self.species_id_ls)} names and ids."
        )

        # filter
        term_list = []
        for species_name in self.species_name_ls:
            term_list.append(f"{species_name}[SCIN]")

        for species_id in self.species_id_ls:
            term_list.append(f"{species_id}[UID]")

        # Get ID
        idlist = []
        for i in range(0, len(term_list), 20):
            retry_times = 0
            while retry_times < 3:
                logger.info(
                    f"Check: {i+1: 4}-{min(i+20, len(self.species_name_ls) + len(self.species_id_ls)): 4}"
                )
                stream = BioEntrez.esearch(
                    db="taxonomy", term=" OR ".join(term_list[i : i + 20])
                )
                record = BioEntrez.read(stream)
                if not record:
                    # not return valid result
                    logger.warning(f"Connect failed. retry times: {retry_times}/3")
                    retry_times += 1
                    continue

                idlist.extend(record["IdList"])
                break

        if not idlist:
            logger.error("Not found any valid species UID")
            return

        logger.info(f"Valid names or ids: {len(idlist)}")

        # Get lineage
        for i in range(0, len(idlist), 20):
            retry_times = 0
            while retry_times < 3:
                logger.info(
                    f"Fetch taxonomy of species: {i+1: 4}-{min(i+20, len(self.species_name_ls) + len(self.species_id_ls)): 4}"
                )
                stream = BioEntrez.efetch(
                    db="taxonomy",
                    id=",".join(idlist[i : i + 20]),
                )
                lineage_records = BioEntrez.read(stream)

                if not lineage_records:
                    # not return valid result
                    logger.warning(f"Connect failed. retry times: {retry_times}/3")
                    retry_times += 1
                    continue

                for r in lineage_records:
                    if r["Rank"] != "species":
                        logger.warning(
                            f"'{r["ScientificName"]}'({r["TaxId"]}) is a {r["Rank"]}, not a species"
                        )

                    lineage_dict = {}

                    for rank in r["LineageEx"]:
                        if rank["Rank"] in ["genus", "family", "order", "class"]:
                            lineage_dict[rank["Rank"]] = rank.get("ScientificName", "?")

                    lineage_dict["species"] = r["ScientificName"]

                    self.lineage_list.append(lineage_dict)
                break

        self._isok = True

    def to_tsv(self, path: str | Path):
        if not self.lineage_list:
            self.get()

        if self._isok:
            df = pd.DataFrame(data=self.lineage_list)
            df.to_csv(path, sep="\t", index=False)


@app_entrez.command(name="species-lineage")
def cmd_species_lineage(
    species_name: str = typer.Option("", "-n", "--species-name"),
    species_id: int = typer.Option(0, "-i", "--species-id"),
):
    if species_name:
        t = GVEUtilsTaxonomy(species_name=species_name)
        t.show()
    elif species_id:
        t = GVEUtilsTaxonomy(species_id=species_id)
        t.show()
    else:
        raise AllFieldsEmptyError("species_name", "species_id")


@app_entrez.command(name="batch-species-lineage")
def cmd_batch_species_lineage(name_or_id_path: str, res_output_path: str):
    species_name_ls, species_id_ls = [], []
    with open(name_or_id_path) as f:
        for line in f:
            line = line.strip()
            try:
                int(line)
            except ValueError:
                species_name_ls.append(line)
                continue
            species_id_ls.append(line)
    t = GVEUtilsBatchTaxonomy(species_name_ls, species_id_ls)
    logger.info("Start searching.")
    t.get()
    t.to_tsv(res_output_path)
