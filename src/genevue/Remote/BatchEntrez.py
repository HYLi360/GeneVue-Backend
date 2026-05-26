from typing import Optional

from Bio import Entrez
from rich.text import Text

from genevue import (
    console,
    AllFieldsEmptyError,
    NothingFoundError,
    ResultIsNotSpeciesError,
)
from genevue.utils.parse import species_italic_name


class XMLParser:
    pass


class NCBIEUtils:
    pass


class ESearch:
    pass


class GeneVueEUtilsTaxonomy:
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

        self.lineage: dict[str, str] = {}

    def get(self):
        subterm_list = []
        for subterm_content, subterm_flag in [
            (self.species_id, "[UID]"),
            (self.species_name, "[SCIN]"),
        ]:
            if subterm_content:
                subterm_list.append("".join([subterm_content, subterm_flag]))

        del subterm_content, subterm_flag

        if not subterm_list:
            console.exception(
                AllFieldsEmptyError(["species_id", "species_name"]),
            )

        console.info(f"Searching by Biopython Entrez Module.")

        # Get ID
        stream = Entrez.esearch(db="taxonomy", term=" AND ".join(subterm_list))
        record = Entrez.read(stream)
        if not record["IdList"]:
            # try to get the true spell
            stream = Entrez.espell(term=self.species_name)
            correct_spell = Entrez.read(stream)["CorrectedQuery"]
            if not self.species_name or not correct_spell:
                console.exception(
                    NothingFoundError(
                        "I tried searching but found nothing. Maybe you should check your input?"
                    )
                )
            else:
                console.exception(
                    NothingFoundError(
                        "I tried searching but found nothing. Maybe you spelled wrong, "
                        f"mixed [{self.species_name}] with [{correct_spell}]?"
                    )
                )

        stream = Entrez.efetch(db="taxonomy", id=record["IdList"][0])
        record = Entrez.read(stream)[0]
        self.species_id = record["TaxId"]
        self.species_name = record["ScientificName"]
        self.genetic_code = record["GeneticCode"]["GCId"]
        self.genetic_code_mt = record["MitoGeneticCode"]["MGCId"]

        if record["Rank"] != "species":
            console.exception(ResultIsNotSpeciesError(record["ScientificName"]))

        for rank in record["LineageEx"]:
            match rank["Rank"]:
                case "genus":
                    self.lineage["genus"] = rank["ScientificName"]
                case "family":
                    self.lineage["family"] = rank["ScientificName"]
                case "order":
                    self.lineage["order"] = rank["ScientificName"]
                case "class":
                    self.lineage["class"] = rank["ScientificName"]

    def show(self):
        console.print(species_italic_name(self.species_name), f"[{self.species_id}]")
        console.print(f"Class:  {self.lineage['class']}")
        console.print(f"Order:  {self.lineage['order']}")
        console.print(f"Family: {self.lineage['family']}")
        console.print(Text("Genus: "), Text(self.lineage["genus"], style="italic"))
        console.print(f"Nucl-Genetic Code Type: {self.genetic_code}")
        console.print(f"MT - Genetic Code Type: {self.genetic_code_mt}")
