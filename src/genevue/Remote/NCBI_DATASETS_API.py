import hashlib
import json
import os
import zipfile
from collections import defaultdict
from enum import Enum
from pathlib import Path
from types import NoneType
from typing import List, Literal, Optional

import requests

from genevue import console
from genevue import local_configure
from genevue.utils.network import Downloader

API_VERSION = "v2"
DATASET_API_BASE_URL = f"https://api.ncbi.nlm.nih.gov/datasets/{API_VERSION}"


class Includes(Enum):
    genome_seq = "GENOME_FASTA"
    gtf = "GENOME_GTF"
    gff3 = "GENOME_GFF"
    gbff = "GENOME_GBFF"
    rna = "RNA_FASTA"
    cds = "CDS_FASTA"
    pep = "PROT_FASTA"
    seq_report = "SEQUENCE_REPORT"


class Datasets4Genome:
    def __init__(
        self,
        accessions: str | List[str],
        chromosomes: Optional[str | List[str]] = None,
        include: List[Includes] | Literal["DEFAULT"] = "DEFAULT",
        hydrated: bool = True,
        chunk_size: int = 20,
        target_dir: Path = Path(".").resolve(),
        zip_name: str = "ncbi_dataset",
        generate_symlinks: bool = False,
        force_make_new_dir: bool = False,
    ):
        self.type = "genome"
        console.info(f"Downloading {self.type} by NCBI Datasets API.")
        self.foldpath = Path(target_dir)
        self.zip_name = zip_name
        self.filepath = Path(target_dir) / f"{zip_name}.zip"

        if force_make_new_dir:
            pathls = list(self.filepath.parents)
            pathls.reverse()
            for path in pathls:
                path.mkdir(exist_ok=True)

        if isinstance(accessions, str):
            self.accession: str = accessions
            self.accession_list: List[str] = [accessions]
        elif isinstance(accessions, list):
            self.accession: str = ", ".join(accessions)
            self.accession_list: List[str] = accessions
        else:
            console.error(
                f"Not supported accession argument type: except 'str' or 'List[str]', but received {type(accessions)}"
            )
        console.info(f"Accession:   {self.accession}")

        if isinstance(chromosomes, str):
            self.chromosome: str = chromosomes
            self.chromosome_list: List[str] = [chromosomes]
        elif isinstance(chromosomes, list):
            self.chromosome: str = ", ".join(chromosomes)
            self.chromosome_list: List[str] = chromosomes
        elif isinstance(chromosomes, NoneType):
            self.chromosome: str = ""
            self.chromosome_list: List[str] = []
        else:
            console.error(
                f"Not supported chromosome argument type: except 'str' or 'List[str]', but received {type(chromosomes)}"
            )
        console.info(
            f"chromosome:  {'Not Specified' if self.chromosome == '' else self.chromosome}"
        )

        if isinstance(include, str):
            self.include: str = include
            self.include_list: List[str] = []
        elif isinstance(include, list):
            self.include: str = ", ".join(
                [include_item.value for include_item in include]
            )
            self.include_list: List[str] = list(
                [include_item.value for include_item in include]
            )
        else:
            console.error(
                f"Not supported include argument type: except 'str' or 'List[str]', but received {type(include)}"
            )
        console.info(f"Include:     {self.include}")

        self.hydrated = "FULLY_HYDRATED" if hydrated else "DATA_REPORT_ONLY"
        console.info(f"Hydrated:    {self.hydrated}")

        self.apikey = local_configure.get_apikey("NCBI")

        if self.apikey is None:
            console.warn("You didn't provide NCBI Datasets API key!")
            console.warn(
                "This isn't a big deal, it only just reduce your maximum RPS from 10 to 5. "
                "However, we still recommend that you apply for an API key."
            )
            console.warn(
                "Visit https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/ for more details."
            )

        if chunk_size > 100:
            raise Exception
        self.chunk_size = chunk_size

        self.generate_symlinks = generate_symlinks

    @property
    def request_download_url(self):
        accessions: str = "%2C".join(self.accession_list)
        chromosomes: str = "&".join(
            [f"chromosomes={chromosome}" for chromosome in self.chromosome_list]
        )
        includes: str = "&".join(
            [f"include_annotation_type={include}" for include in self.include_list]
        )
        hydrated: str = (
            "" if self.hydrated != "DATA_REPORT_ONLY" else f"hydrated={self.hydrated}"
        )
        filename: str = (
            ""
            if self.filepath.name == "ncbi_dataset.zip"
            else f"filename={self.zip_name}.zip"
        )
        params = "&".join(
            [item for item in [chromosomes, includes, hydrated, filename] if item != ""]
        )

        return f"{DATASET_API_BASE_URL}/{self.type}/accession/{accessions}/download?{params}"

    @property
    def request_check_url(self):
        accessions: str = "%2C".join(self.accession_list)
        return f"{DATASET_API_BASE_URL}/{self.type}/accession/{accessions}/check"

    @property
    def request_header(self):
        header = {"accept": "application/json"}
        if self.apikey is not None:
            header["api_key"] = self.apikey

        return header

    @staticmethod
    def _chunk_iter(data, chunk_size):
        for i in range(0, len(data), chunk_size):
            yield data[i : i + chunk_size]

    def _check_invalid_duplicated_accessions(self) -> None:
        console.info("Start Check.")
        console.info(f"Input accession: {len(self.accession_list)}")
        need_update_accession_list_label = False

        # Check Duplicated
        new_accession_list = list(dict.fromkeys(self.accession_list))
        console.info(
            f"Accessions count after remove duplicated: {len(new_accession_list)}"
        )
        if len(new_accession_list) != len(self.accession_list):
            console.warn("Duplicate accession detected.")
            need_update_accession_list_label = True

        # Check Invalid
        console.info(f"Try opening the url {self.request_check_url}")
        valid_accessions, invalid_accessions = [], []

        for chunk in self._chunk_iter(new_accession_list, 100):
            self.accession_list = chunk
            res = requests.get(
                self.request_check_url, headers=self.request_header
            ).json()
            valid_accessions.extend(res.get("valid_assemblies", []))
            invalid_accessions.extend(res.get("invalid_assemblies", []))

        if len(valid_accessions) != len(new_accession_list):
            console.warn("Invalid accession detected.")
            console.warn(f"Invalid accessions: {invalid_accessions}")
            need_update_accession_list_label = True
        else:
            console.info("Accessions check passed.")

        if need_update_accession_list_label:
            console.warn("Accession list updated.")
        self.accession_list = valid_accessions

    def download(self):
        # Download
        accession_list_original = self.accession_list

        md5file_path = self.foldpath / "md5sum.txt"
        asm_data_report_path = (
            self.foldpath / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
        )
        dataset_catalog_path = (
            self.foldpath / "ncbi_dataset" / "data" / "dataset_catalog.json"
        )
        dataset_catalog_list = []
        md5dict = defaultdict(str)
        for order, chunk in enumerate(
            self._chunk_iter(accession_list_original, self.chunk_size)
        ):
            console.info(f"Start download batch {order+1}.")
            self.accession_list = chunk

            # check
            self._check_invalid_duplicated_accessions()

            # download
            downloader = Downloader(
                self.request_download_url,
                self.request_header,
                self.foldpath
                / f"{self.zip_name}.{order+1}.zip",  # avoid it to be covered
            )
            downloader.download()

            # extract.
            with zipfile.ZipFile(
                self.foldpath / f"{self.zip_name}.{order+1}.zip", "r"
            ) as zipf:
                console.info("Extracting.")
                zipf.extractall(path=self.foldpath)

            # read md5.
            with open(md5file_path) as md5f:
                for line in md5f:
                    md5sum, filepath = line.strip().split()
                    if filepath == "ncbi_dataset/data/assembly_data_report.jsonl":
                        md5dict[
                            f"ncbi_dataset/data/assembly_data_report.{order+1}.jsonl"
                        ] = md5sum
                    elif filepath == "ncbi_dataset/data/dataset_catalog.json":
                        md5dict[f"ncbi_dataset/data/dataset_catalog.{order+1}.json"] = (
                            md5sum
                        )
                    else:
                        md5dict[filepath] = md5sum

            # concat jsonl.
            with (
                open(asm_data_report_path) as asm_original,
                open(
                    self.foldpath
                    / "ncbi_dataset"
                    / "data"
                    / "assembly_data_report_concat.jsonl",
                    "a",
                ) as asm_final,
            ):
                for line in asm_original:
                    asm_final.write(line)

            # read the dataset catalog list.
            # original format:
            # position 0: assembly_data_report.jsonl (no needed)
            # position 1: accession X
            # position 2: accession Y
            # ...

            # {
            # "apiVersion": "V2",
            # "assemblies": [
            # {
            #   "files": [
            #     {
            #       "filePath": "assembly_data_report.jsonl",
            #       "fileType": "DATA_REPORT",
            #       "uncompressedLengthBytes": "54774"
            #     }
            #   ]
            # },{
            #   "accession": "GCA_028411795.1",
            #   "files": [
            #     {
            #       "filePath": "GCA_028411795.1/genomic.gtf",
            #       "fileType": "GTF",
            #       "uncompressedLengthBytes": "116499833"
            #     },
            #     {
            #       "filePath": "GCA_028411795.1/protein.faa",
            #       "fileType": "PROTEIN_FASTA",
            #       "uncompressedLengthBytes": "13042297"
            #     }
            #   ]
            # }, ...
            # ]}
            dataset_catalog_list.extend(
                json.load(
                    open(
                        self.foldpath / "ncbi_dataset" / "data" / "dataset_catalog.json"
                    )
                )["assemblies"][1:]
            )

            # rename.
            os.rename(md5file_path, self.foldpath / f"md5sum.{order+1}.txt")
            os.rename(
                asm_data_report_path,
                self.foldpath
                / "ncbi_dataset"
                / "data"
                / f"assembly_data_report.{order+1}.jsonl",
            )
            os.rename(
                dataset_catalog_path,
                self.foldpath
                / "ncbi_dataset"
                / "data"
                / f"dataset_catalog.{order+1}.json",
            )

        # write back the dataset_catalog.json and md5 dict
        # just need some rebuild to make it much easier to read...

        # Unified the file types.
        # possible file types:
        # - "CDS_NUCLEOTIDE_FASTA" -> cds
        # - "GENBANK_FLAT_FILE" -> gbff
        # - "GENOMIC_NUCLEOTIDE_FASTA" -> genome_seq
        # - "GFF3" -> "gff3"
        # - "GTF" -> "gtf"
        # - "PROTEIN_FASTA" -> pep
        # - "RNA_NUCLEOTIDE_FASTA" -> rna
        # - "SEQUENCE_REPORT" -> seq_report
        type_switch_dict_type = {
            "CDS_NUCLEOTIDE_FASTA": "CDS_FASTA",
            "GENBANK_FLAT_FILE": "GENOME_GBFF",
            "GENOMIC_NUCLEOTIDE_FASTA": "GENOME_FASTA",
            "GFF3": "GENOME_GTF",
            "GTF": "GENOME_GFF",
            "PROTEIN_FASTA": "PROT_FASTA",
            "RNA_NUCLEOTIDE_FASTA": "RNA_FASTA",
            "SEQUENCE_REPORT": "SEQUENCE_REPORT",
        }

        # Flatten the path index
        # to
        # {accession: {file_type: file_path}}
        dataset_catalog_simplified = {}
        for accession_block in dataset_catalog_list:
            dataset_catalog_simplified[accession_block["accession"]] = {}
            for file in accession_block["files"]:
                dataset_catalog_simplified[accession_block["accession"]][
                    type_switch_dict_type[file["fileType"]]
                ] = file["filePath"]
        json.dump(
            dataset_catalog_simplified,
            open(
                self.foldpath
                / "ncbi_dataset"
                / "data"
                / f"dataset_catalog_simplified.json",
                "w",
            ),
        )
        # for md5...
        with open(self.foldpath / "md5sum.txt", "w") as md5f:
            for path in md5dict:
                md5f.write(f"{path} {md5dict[path]}\n")

        # restore self.accession_list
        self.accession_list = accession_list_original

        # check md5 and generate symlinks (if needed.)
        for accession in dataset_catalog_simplified:
            files_dict = dataset_catalog_simplified[accession]
            for include in self.include_list:
                try:
                    sub_file_path = (
                        self.foldpath / "ncbi_dataset" / "data" / files_dict[include]
                    )
                except KeyError:
                    console.warn(f"Accession {accession} has no file of {include}")
                    continue

                # check md5.
                with open(sub_file_path, "rb") as f:
                    console.info(f"Checking {sub_file_path}")
                    md5sum_hex = hashlib.file_digest(f, "md5").hexdigest()
                    if (
                        md5dict[f"ncbi_dataset/data/{files_dict[include]}"]
                        != md5sum_hex
                    ):
                        console.warn(
                            f"Wrong md5 sum of {sub_file_path}! except {md5dict[files_dict[include]]}, get {md5sum_hex}"
                        )
                        continue

                # create symlink.
                if self.generate_symlinks:
                    console.info(f"Set symlink for {sub_file_path}")
                    target = self.foldpath / "symlinks" / include / accession
                    target.parent.parent.mkdir(exist_ok=True)
                    target.parent.mkdir(exist_ok=True)
                    os.symlink(src=sub_file_path, dst=target)
