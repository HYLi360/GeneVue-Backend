import asyncio
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

from genevue import setup_rich_logger, console
from genevue.utils.network import AsyncDownloadManager

logger = setup_rich_logger(__name__, console)

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
        apikey: Optional[str] = None,
        generate_symlinks: bool = False,
        force_make_new_dir: bool = False,
        max_concurrent: int = 3,
    ):
        self.type = "genome"
        logger.info(f"Downloading {self.type} by NCBI Datasets API.")
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
            logger.error(
                f"Not supported accession argument type: except 'str' or 'List[str]', but received {type(accessions)}"
            )
        logger.info(f"Accession:   {self.accession}")

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
            logger.error(
                f"Not supported chromosome argument type: except 'str' or 'List[str]', but received {type(chromosomes)}"
            )
        logger.info(
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
            logger.error(
                f"Not supported include argument type: except 'str' or 'List[str]', but received {type(include)}"
            )
        logger.info(f"Include:     {self.include}")

        self.hydrated = "FULLY_HYDRATED" if hydrated else "DATA_REPORT_ONLY"
        logger.info(f"Hydrated:    {self.hydrated}")

        self.apikey = apikey
        if self.apikey is None:
            logger.warning(
                "You didn't provide NCBI Datasets API key!\n"
                "This isn't a big deal, it only just reduce your maximum RPS from 10 to 5. "
                "However, we still recommend that you apply for an API key.\n"
                "Visit https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/ for more details."
            )
            if max_concurrent < 0 or max_concurrent > 5:
                logger.warning(
                    "Illegal value for max_concurrent argument! Set to 3 automatically."
                )
                self.max_concurrent = 3
            else:
                self.max_concurrent = max_concurrent
        else:
            if max_concurrent < 0 or max_concurrent > 10:
                logger.warning(
                    "Illegal value for max_concurrent argument! Set to 3 automatically."
                )
                self.max_concurrent = 3
            else:
                self.max_concurrent = max_concurrent

        if chunk_size > 100:
            logger.warning("Too much chunk size! Automately set to 20.")
        self.chunk_size = chunk_size

        self.generate_symlinks = generate_symlinks

    # URL builders
    def build_download_url(self, accession_list: List[str]) -> str:
        accessions: str = "%2C".join(accession_list)
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

    def build_check_url(self, accession_list: List[str]) -> str:
        accessions: str = "%2C".join(accession_list)
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

    # Validation
    def _check_invalid_duplicated_accessions(
        self, accession_list: List[str]
    ) -> List[str]:
        logger.info("Start check.")
        logger.info(f"Input accession: {len(accession_list)}")

        new_list = list(dict.fromkeys(accession_list))
        logger.info(f"Accessions count after remove duplicated: {len(new_list)}")
        if len(new_list) != len(accession_list):
            logger.warning("Duplicate accession detected.")

        logger.info("Validating accessions via NCBI check API...")
        valid_accessions, invalid_accessions = [], []

        for chunk in self._chunk_iter(new_list, 100):
            url = self.build_check_url(chunk)
            logger.info(f"Try opening the url {url}")
            res = requests.get(url, headers=self.request_header).json()
            valid_accessions.extend(res.get("valid_assemblies", []))
            invalid_accessions.extend(res.get("invalid_assemblies", []))

        if len(valid_accessions) != len(new_list):
            logger.warning("Invalid accession detected.")
            logger.warning(f"Invalid accessions: {invalid_accessions}")
        else:
            logger.info("Accessions check passed.")
        return valid_accessions

    def _validate_and_chunk(self, accessions: List[str]) -> List[List[str]]:
        valid = self._check_invalid_duplicated_accessions(accessions)
        if not valid:
            return []
        return list(self._chunk_iter(valid, self.chunk_size))

    # Post-download processing
    def _process_downloaded_chunk(
        self, order: int, dataset_catalog_list: list, md5dict: dict
    ):
        logger.info(f"Processing batch {order+1}.")

        md5file_path = self.foldpath / "md5sum.txt"
        asm_data_report_path = (
            self.foldpath / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
        )
        dataset_catalog_path = (
            self.foldpath / "ncbi_dataset" / "data" / "dataset_catalog.json"
        )

        # extract
        with zipfile.ZipFile(
            self.foldpath / f"{self.zip_name}.{order+1}.zip", "r"
        ) as zipf:
            logger.info("Extracting.")
            zipf.extractall(path=self.foldpath)

        # read md5
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

        # concat jsonl
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

        # read the dataset catalog list
        dataset_catalog_list.extend(
            json.load(open(dataset_catalog_path))["assemblies"][1:]
        )

        # rename
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
            self.foldpath / "ncbi_dataset" / "data" / f"dataset_catalog.{order+1}.json",
        )

    # finalization
    def _finalize_download(self, dataset_catalog_list: list, md5dict: dict):
        type_switch_dict_type = {
            "CDS_NUCLEOTIDE_FASTA": "CDS_FASTA",
            "GENBANK_FLAT_FILE": "GENOME_GBFF",
            "GENOMIC_NUCLEOTIDE_FASTA": "GENOME_FASTA",
            "GFF3": "GENOME_GFF",
            "GTF": "GENOME_GTF",
            "PROTEIN_FASTA": "PROT_FASTA",
            "RNA_NUCLEOTIDE_FASTA": "RNA_FASTA",
            "SEQUENCE_REPORT": "SEQUENCE_REPORT",
        }

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
                / "dataset_catalog_simplified.json",
                "w",
            ),
        )

        # write back the md5 dict
        with open(self.foldpath / "md5sum.txt", "w") as md5f:
            for path in md5dict:
                md5f.write(f"{path} {md5dict[path]}\n")

        # check md5 and generate symlinks (if needed)
        for accession in dataset_catalog_simplified:
            files_dict = dataset_catalog_simplified[accession]
            for include in self.include_list:
                try:
                    sub_file_path = (
                        self.foldpath / "ncbi_dataset" / "data" / files_dict[include]
                    )
                except KeyError:
                    logger.warning(f"Accession {accession} has no file of {include}")
                    continue

                with open(sub_file_path, "rb") as f:
                    logger.info(f"Checking {sub_file_path}")
                    md5sum_hex = hashlib.file_digest(f, "md5").hexdigest()
                    if (
                        md5dict[f"ncbi_dataset/data/{files_dict[include]}"]
                        != md5sum_hex
                    ):
                        logger.warning(
                            f"Wrong md5 sum of {sub_file_path}! "
                            f"except {md5dict[files_dict[include]]}, get {md5sum_hex}"
                        )
                        continue

                if self.generate_symlinks:
                    logger.info(f"Set symlink for {sub_file_path}")
                    target = self.foldpath / "symlinks" / include / accession
                    target.parent.parent.mkdir(exist_ok=True)
                    target.parent.mkdir(exist_ok=True)
                    os.symlink(src=sub_file_path, dst=target)

    # Parallel download orchestration
    async def _download_chunks_parallel(
        self, chunks: List[List[str]], dataset_catalog_list: list, md5dict: dict
    ):
        semaphore = asyncio.Semaphore(self.max_concurrent)

        async def download_one(order: int, chunk: List[str]):
            url = self.build_download_url(chunk)
            zip_path = self.foldpath / f"{self.zip_name}.{order+1}.zip"
            logger.info(f"Start download batch {order+1}.")
            async with semaphore:
                dm = AsyncDownloadManager(url, self.request_header, zip_path)
                await dm.download()

        tasks = [download_one(i, c) for i, c in enumerate(chunks)]
        await asyncio.gather(*tasks)

        for order in range(len(chunks)):
            self._process_downloaded_chunk(order, dataset_catalog_list, md5dict)

    # Public API
    def download(self):
        accession_list_original = self.accession_list

        chunks = self._validate_and_chunk(accession_list_original)
        if not chunks:
            logger.error("No valid accessions to download.")
            return

        dataset_catalog_list = []
        md5dict = defaultdict(str)

        if len(chunks) == 1:
            url = self.build_download_url(chunks[0])
            logger.info("Start download batch 1.")
            dm = AsyncDownloadManager(
                url,
                self.request_header,
                self.foldpath / f"{self.zip_name}.1.zip",
            )
            asyncio.run(dm.download())
            self._process_downloaded_chunk(0, dataset_catalog_list, md5dict)
        else:
            logger.info(
                f"Downloading {len(chunks)} batches in parallel "
                f"(max {self.max_concurrent} concurrent)."
            )
            asyncio.run(
                self._download_chunks_parallel(chunks, dataset_catalog_list, md5dict)
            )

        self.accession_list = accession_list_original
        self._finalize_download(dataset_catalog_list, md5dict)
