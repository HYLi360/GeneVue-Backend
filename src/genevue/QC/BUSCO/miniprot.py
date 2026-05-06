# coding: utf-8
"""
miniprot.py

Module for running Miniprot.

Author(s): Matthew Berkeley

Copyright (c) 2015-2025, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from genevue.QC.BUSCO.base import GenePredictor, BaseRunner
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
import subprocess
import os
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import re
from multiprocessing import Pool

logger = BuscoLogger.get_logger(__name__)


class MiniprotRunner(BaseRunner):
    """
    Class to run Miniprot
    """

    name = "miniprot"
    cmd = "miniprot"

    def __init__(self):
        super().__init__()
        self._output_folder = os.path.join(self.run_folder, "miniprot_output")
        self.translated_proteins_folder = os.path.join(
            self._output_folder, "translated_proteins"
        )
        self.create_dirs(
            [
                self._output_folder,
                self.translated_proteins_folder,
            ]
        )
        self.index_file = os.path.join(self._output_folder, "ref.mpi")
        self.refseq_db = None
        self.incomplete_buscos = None
        self._output_basename = None

    def configure_runner(self, *args):
        """
        Configure Miniprot runner
        """
        super().configure_runner(*args)
        self.run_number += 1

    def check_tool_dependencies(self):
        pass

    def configure_job(self, *args):
        """
        Overridden by child classes
        """
        return

    def generate_job_args(self):
        yield

    def get_version(self):
        help_output = subprocess.check_output(
            [self.cmd, "--version"], stderr=subprocess.STDOUT, shell=False
        )
        version = help_output.decode("utf-8").strip()
        return version

    @property
    def output_folder(self):
        return self._output_folder

    def reset(self):
        super().reset()

    def run(self):
        super().run()
        self.total = 1
        self.run_jobs()


class MiniprotIndexRunner(MiniprotRunner):
    name = "miniprot_index"

    def generate_job_args(self):
        yield "index"

    def configure_job(self, *args):
        """
        Configure Miniprot job
        """

        miniprot_job = self.create_job()
        miniprot_job.add_parameter("-t")
        miniprot_job.add_parameter(str(self.cpus))
        miniprot_job.add_parameter("-d")
        miniprot_job.add_parameter(self.index_file)
        miniprot_job.add_parameter(self.input_file)

        return miniprot_job


class MiniprotAlignRunner(MiniprotRunner, GenePredictor):
    name = "miniprot_align"

    def __init__(self):
        super().__init__()
        self.output_gff = None

        self.sequences_aa = {}
        self.busco_matches = defaultdict(set)
        self.gene_matches = defaultdict(list)
        self.combined_pred_protein_seqs = os.path.join(
            self._output_folder, "combined_pred_proteins.fas"
        )
        self.output_sequences = []
        self.gene_nominal = 0
        self.gene_lookup = {}
        self.cigar_lookup = {}
        self.nominal_lookup = defaultdict(list)

        self.gene_details = defaultdict(dict)
        self.sequences_aa = {}
        self.sequences_nt = {}

    def generate_job_args(self):
        yield "align"

    def configure_job(self, *args):

        miniprot_job = self.create_job()
        miniprot_job.add_parameter("--trans")
        miniprot_job.add_parameter("-u")
        miniprot_job.add_parameter("-I")
        miniprot_job.add_parameter("--outs")
        miniprot_job.add_parameter("0.95")
        miniprot_job.add_parameter("-t")
        miniprot_job.add_parameter(str(self.cpus))
        miniprot_job.add_parameter("--gff")

        miniprot_job.add_parameter(self.index_file)
        miniprot_job.add_parameter(self.refseq_db)

        return miniprot_job

    def configure_runner(self, incomplete_buscos=None):
        super().configure_runner([])
        self.logfile_path_out = os.path.join(
            self.config.get("busco_run", "main_out"),
            "logs",
            "{}_{}_out.log".format(self.name, os.path.basename(self.lineage_dataset)),
        )
        self.logfile_path_err = (
            self.logfile_path_out.rpartition("_out.log")[0] + "_err.log"
        )

        self.incomplete_buscos = incomplete_buscos
        self._output_basename = os.path.join(
            self._output_folder, os.path.basename(self.input_file)
        )
        gzip_refseq = os.path.join(self.lineage_dataset, "refseq_db.faa.gz")
        self.refseq_db = self.decompress_refseq_file(gzip_refseq)
        self.output_gff = Path(self._output_folder).joinpath(
            "{}_{}{}".format(
                Path(self.input_file).stem,
                os.path.basename(self.lineage_dataset),
                ".gff",
            )
        )

    def create_symlink(self):
        if not self.output_gff.exists():
            Path(self.output_gff).symlink_to(
                os.path.relpath(self.logfile_path_out, self.output_gff.parent)
            )
        return

    def save_record(
        self,
        gene_id,
        target_id,
        contig_id,
        contig_start,
        contig_end,
        strand,
        score,
        exon_coords,
        ata_seq,
        protein_start,
        protein_end,
        protein_length,
        stop_codon_count,
        identity,
        gff_start,
        gff_end,
    ):
        self.all_gff_records.append(
            np.array(
                (
                    gene_id,
                    target_id.split("_")[0],
                    target_id,
                    contig_id,
                    contig_start,
                    contig_end,
                    strand,
                    score,
                    exon_coords,
                    ata_seq,
                    protein_start,
                    protein_end,
                    protein_length,
                    stop_codon_count,
                    identity,
                    gff_start,
                    gff_end,
                ),
                dtype=[
                    ("gene_id", "U500"),
                    ("busco_id", "U100"),
                    ("target_id", "U500"),
                    ("contig_id", "U500"),
                    ("contig_start", "i4"),
                    ("contig_end", "i4"),
                    ("strand", "U1"),
                    ("score", "i4"),
                    ("exon_coords", "O"),
                    ("aa_seq", "U10000"),
                    ("protein_start", "i4"),
                    ("protein_end", "i4"),
                    ("protein_length", "i4"),
                    ("stop_codon_count", "i4"),
                    ("identity", "f4"),
                    ("gff_start", "i8"),
                    ("gff_end", "i8"),
                ],
            )
        )
        return

    def parse_output(self):
        self.create_symlink()
        self.all_gff_records = []
        hits_by_group = defaultdict(list)
        with open(self.output_gff, "r") as gff:
            line = gff.readline()
            # Load all gff record headers
            while line:
                if line.startswith("##PAF"):
                    cursor_position = gff.tell() - len(line)
                    fields = line.strip().split("\t")[1:]
                    if fields[5] == "*":
                        ## Unmapped protein
                        line = gff.readline()
                        continue

                    target_id = fields[0]
                    contig_id = fields[5]
                    contig_start = int(fields[7])
                    contig_end = int(fields[8])

                    score = int(fields[13].strip().split(":")[2])

                    hits_by_group[target_id.split("_")[0]].append(
                        (
                            target_id,
                            score,
                            contig_id,
                            contig_start,
                            contig_end,
                            cursor_position,
                        )
                    )

                line = gff.readline()

            # Filter gff records, removing any with greater than 80% overlap with highest scoring region
            all_regions = []
            for group, hits in hits_by_group.items():
                hits.sort(key=lambda x: -x[1])
                selected_regions = []
                for hit in hits:
                    (
                        target_id,
                        score,
                        contig_id,
                        contig_start,
                        contig_end,
                        cursor_position,
                    ) = hit
                    l1 = contig_end - contig_start
                    for entry in selected_regions:
                        _, _, c_id, s, e, _ = entry
                        if c_id == contig_id and (
                            (contig_start <= s <= contig_end)
                            or (s <= contig_start <= e)
                        ):
                            l2 = e - s
                            l_overlap = min(e, contig_end) - max(s, contig_start)
                            if l_overlap > 0.8 * max(l1, l2):
                                break

                    else:
                        selected_regions.append(
                            (
                                target_id,
                                score,
                                contig_id,
                                contig_start,
                                contig_end,
                                cursor_position,
                            )
                        )

                all_regions.extend(selected_regions)

            # Parse only the records that passed the filter
            gff.seek(0)
            for region in all_regions:
                (
                    target_id,
                    score,
                    contig_id,
                    contig_start,
                    contig_end,
                    cursor_position,
                ) = region
                gff.seek(cursor_position)
                gff_start = gff_end = None
                line = gff.readline()
                fields = line.strip().split("\t")[1:]
                protein_length = int(fields[1])
                protein_start = int(fields[2])
                protein_end = int(fields[3])
                strand = fields[4]
                gene_id = "{}|{}:{}-{}|{}".format(
                    target_id, contig_id, contig_start, contig_end, strand
                )
                stop_codon_count = int(fields[17].strip().split(":")[2])
                exon_coords = []

                sta_line = gff.readline()
                sta_seq = sta_line.strip().split("\t")[1]
                ata_seq = sta_seq.upper()
                mrna_line = gff.readline()
                if "mRNA" in mrna_line:
                    gff_start = gff.tell() - len(mrna_line)
                    info_str = mrna_line.split("\t")[8]
                    info_dict = dict(v.split("=") for v in info_str.split(";"))
                    identity = float(info_dict["Identity"])
                    exons_line = gff.readline()
                    while exons_line and not exons_line.startswith("##PAF"):
                        fields = exons_line.strip().split("\t")
                        if fields[2] == "CDS":
                            start, stop, score, strand = (
                                int(fields[3]),
                                int(fields[4]),
                                float(fields[5]),
                                fields[6],
                            )
                            exon_coords.append((start, stop, score, strand))
                        exons_line = gff.readline()
                    gff_end = gff.tell() - len(exons_line)

                else:
                    identity = 0.0

                self.save_record(
                    gene_id,
                    target_id,
                    contig_id,
                    contig_start,
                    contig_end,
                    strand,
                    score,
                    exon_coords,
                    ata_seq,
                    protein_start,
                    protein_end,
                    protein_length,
                    stop_codon_count,
                    identity,
                    gff_start,
                    gff_end,
                )

        self.filtered_matches = np.array(self.all_gff_records)
        return

    def record_gene_details(self):
        for match in self.filtered_matches:
            gene_id = match["gene_id"].item()
            target_id = match["target_id"].item()
            contig_id = match["contig_id"].item()
            contig_start = match["contig_start"].item()
            contig_end = match["contig_end"].item()
            strand = match["strand"].item()
            score = match["score"].item()
            exon_coords = match["exon_coords"]
            ata_seq = match["aa_seq"].item()
            protein_start = match["protein_start"].item()
            protein_end = match["protein_end"].item()
            protein_length = match["protein_length"].item()
            stop_codon_count = match["stop_codon_count"].item()
            identity = match["identity"].item()
            gff_start = match["gff_start"].item()
            gff_end = match["gff_end"].item()

            self.gene_details[gene_id].update(
                {
                    "gene_id": gene_id,
                    "target_id": target_id,
                    "contig_id": contig_id,
                    "contig_start": contig_start,
                    "contig_end": contig_end,
                    "strand": strand,
                    "score": score,
                    "exon_coords": exon_coords,
                    "aa_seq": SeqRecord(Seq(ata_seq), id=gene_id),
                    "protein_start": protein_start,
                    "protein_end": protein_end,
                    "protein_length": protein_length,
                    "stop_codon_count": stop_codon_count,
                    "identity": identity,
                    "run_number": self.run_number,
                    "gff_start": gff_start,
                    "gff_end": gff_end,
                }
            )
            self.sequences_aa[gene_id] = SeqRecord(Seq(ata_seq), id=gene_id)
            self.busco_matches[target_id.split("_")[0]].add(gene_id)
        return

    def check_overlap(self, contig_id, contig_start, contig_end, score, gene_id):
        keeper = gene_id
        matches_to_remove = []
        matches = np.array(self.genomic_regions[contig_id])
        overlap_matches = matches[
            (
                (matches["contig_start"] >= contig_start)
                & (matches["contig_start"] < contig_end)
            )
            | (
                (matches["contig_end"] > contig_start)
                & (matches["contig_end"] <= contig_end)
            )
        ]
        if len(overlap_matches) > 0:
            for match in overlap_matches:
                overlap_start = max(contig_start, match["contig_start"])
                overlap_end = min(contig_end, match["contig_end"])
                overlap_length = overlap_end - overlap_start
                if (overlap_length > 0.8 * (contig_end - contig_start)) or (
                    overlap_length > 0.8 * (match["contig_end"] - match["contig_start"])
                ):
                    if score > match["score"]:
                        keeper = gene_id
                        matches_to_remove.append(match)
                    else:
                        keeper = match["gene_id"]
        if len(matches_to_remove) > 0:
            self.genomic_regions[contig_id] = [
                m for m in matches if m not in matches_to_remove
            ]
        return keeper

    def write_protein_sequences_per_busco(self):
        for busco_id in self.busco_matches:
            seqs_to_write = []
            output_filename = os.path.join(
                self.translated_proteins_folder, "{}.faa".format(busco_id)
            )
            self.output_sequences.append(output_filename)
            with open(output_filename, "w") as f:
                for g in self.busco_matches[busco_id]:
                    if g in self.sequences_aa:
                        seqs_to_write.append(self.sequences_aa[g])
                SeqIO.write(seqs_to_write, f, "fasta")
