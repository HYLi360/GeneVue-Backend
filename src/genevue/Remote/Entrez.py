from typing import Optional, Annotated

import typer
from Bio import Entrez as BioEntrez
from rich.text import Text

from genevue import (
    setup_rich_logger,
    console,
    AllFieldsEmptyError,
    NothingFoundError,
    ResultIsNotSpeciesError,
)
from genevue.configure import Configure
from genevue.utils.parse import species_italic_name

logger = setup_rich_logger(__name__, console)

app_entrez = typer.Typer()


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

        configure = Configure()

        BioEntrez.email = configure.email

        if not subterm_list:
            logger.exception(
                AllFieldsEmptyError(["species_id", "species_name"]),
            )

        logger.info(f"Searching by Biopython Entrez Module.")

        # Get ID
        stream = BioEntrez.esearch(db="taxonomy", term=" AND ".join(subterm_list))
        record = BioEntrez.read(stream)
        if not record["IdList"]:
            # try to get the true spell
            stream = BioEntrez.espell(term=self.species_name)
            correct_spell = BioEntrez.read(stream)["CorrectedQuery"]
            if not self.species_name or not correct_spell:
                logger.exception(
                    NothingFoundError(
                        "I tried searching but found nothing. Maybe you should check your input?"
                    )
                )
            else:
                logger.exception(
                    NothingFoundError(
                        "I tried searching but found nothing. Maybe you spelled wrong, "
                        f"mixed [{self.species_name}] with [{correct_spell}]?"
                    )
                )

        stream = BioEntrez.efetch(db="taxonomy", id=record["IdList"][0])
        record = BioEntrez.read(stream)[0]
        self.species_id = record["TaxId"]
        self.species_name = record["ScientificName"]
        self.genetic_code = record["GeneticCode"]["GCId"]
        self.genetic_code_mt = record["MitoGeneticCode"]["MGCId"]

        if record["Rank"] != "species":
            logger.exception(ResultIsNotSpeciesError(record["ScientificName"]))

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


@app_entrez.command(name="species-lineage")
def cmd_species_lineage(
    species_name: Annotated[str, typer.Argument()] = "",
):
    t = GeneVueEUtilsTaxonomy(species_name=species_name)
    t.get()
    t.show()
