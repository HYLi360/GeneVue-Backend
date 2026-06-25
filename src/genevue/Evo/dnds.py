#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES


# ---Import.------------------------------------------------------------------------------------------
from multiprocessing import Pool
from rich.progress import Progress
from typing import Literal
from math import ceil
import time
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from genevue import console
from genevue.External.PAML import YN00, ConfigCODEML, CODEML
from genevue.External.frontend import muscle5, mafft


# ----------------------------------------------------------------------------------------------------
def _codon_align(
    cds_records: tuple[SeqRecord, SeqRecord],
    pep_aligned: tuple[SeqRecord, SeqRecord],
) -> str:
    # init
    cds1, cds2 = str(cds_records[0].seq), str(cds_records[1].seq)
    pep1, pep2 = (
        str(pep_aligned[0].seq),
        str(pep_aligned[1].seq),
    )
    cdsp1, cdsp2 = 0, 0
    codons_buf1, codons_buf2 = [], []

    # len(pep1) == len(pep2)
    for i in range(len(pep1)):
        pep_site1, pep_site2 = pep1[i], pep2[i]
        # Extract two codons when...
        # pep1 ......A......
        # pep2 ......A......
        if (pep_site1 != "-") and (pep_site2 != "-"):
            # pick codons
            codon1, codon2 = (
                cds1[cdsp1 : cdsp1 + 3],
                cds2[cdsp2 : cdsp2 + 3],
            )
            # concat codons
            codons_buf1.append(codon1)
            codons_buf2.append(codon2)
            # move the pointers
            cdsp1 += 3
            cdsp2 += 3
        # Skip when...
        # pep1 ......-......
        # pep2 ......A......
        elif pep_site1 != "-":
            cdsp1 += 3
        elif pep_site2 != "-":
            cdsp2 += 3
        continue

    # generate paml input
    codons_buf1 = "".join(codons_buf1)
    codons_buf2 = "".join(codons_buf2)
    # Structure of paml input: header (sequence count (4-char length) + paired base count (7-char length)
    # seq1
    # (Detailed sequence)
    # seq2
    # (Detailed sequence)
    # Concat with "\n" and join method.
    header = f"{2: 4}{len(codons_buf1): 7}"
    return "\n".join([header, "seq1", codons_buf1, "seq2", codons_buf2])


# ---Dispatch by algorithm/executable program name (MUSCLE/MAFFT).------------------------------------
def _alignp(
    method: Literal["muscle", "mafft", None],
    workdir: str,
    pep_file_name: str,
) -> tuple[SeqRecord, SeqRecord]:
    match method:
        case "muscle":
            muscle5(workdir, pep_file_name, "align.faa")
        case "mafft":
            mafft(workdir, pep_file_name, "align.faa")
    return tuple(SeqIO.parse(f"{workdir}/align.faa", "fasta"))  # type: ignore


# ----------------------------------------------------------------------------------------------------
def _calc_yn00_tasks(
    cds: tuple[SeqRecord, SeqRecord],
    pep: tuple[SeqRecord, SeqRecord],
    align_method: Literal["muscle", "mafft"],
    yn00_p: YN00,
) -> dict:
    # Tasks should focus on executing the task at hand, rather than preparing in advance to handle exceptions.
    # Create a temporary folder for single tasks.
    with tempfile.TemporaryDirectory(dir="/dev/shm") as tmpdir:
        # Step0: Write records as file.
        SeqIO.write(cds, f"{tmpdir}/cds.fna", "fasta")

        # Step1: pep Alignment.
        SeqIO.write(pep, f"{tmpdir}/pep.faa", "fasta")
        pep_aligned = _alignp(align_method, tmpdir, "pep.faa")

        # Step2: Codon Alignment.
        codon_aln = _codon_align(cds, pep_aligned)

        # Step3: do yn00.
        result = yn00_p.run(tmpdir, codon_aln)
        result["seq1"] = cds[0].id
        result["seq2"] = cds[1].id

        return result


def calc_yn00(
    cds1: str | list[SeqRecord],
    cds2: str | list[SeqRecord],
    pep1: str | list[SeqRecord],
    pep2: str | list[SeqRecord],
    align_method: Literal["muscle", "mafft"],
    yn00_p: YN00,
    processes_num: int = 2,
) -> list:
    """Calculate the dN and dS values between homologous gene pairs by yn00."""
    yn00_start = time.time()
    total_records = []
    for records in [cds1, cds2, pep1, pep2]:
        total_records.append(
            [record for record in SeqIO.parse(records, "fasta")]
            if isinstance(records, str)
            else records
        )

    cds, pep = list(zip(total_records[0], total_records[1])), list(
        zip(total_records[2], total_records[3])
    )

    args = [(cds[i], pep[i], align_method, yn00_p) for i in range(len(cds))]
    chunks = max(1, ceil(len(args) / processes_num))

    with Progress() as progress:
        task = progress.add_task("YN00 running.", total=len(args))

        def _cb(_):
            progress.update(task, advance=1, chunks=chunks)

        with Pool(processes=processes_num) as pool:
            async_results = [
                pool.apply_async(_calc_yn00_tasks, a, callback=_cb) for a in args
            ]
            results = [ar.get() for ar in async_results]

    used_seconds = time.time() - yn00_start
    console.log(
        f"YN00 finished in {used_seconds//60:.0f} min {used_seconds%60:.1f} sec."
    )
    return results


# ----------------------------------------------------------------------------------------------------
def _calc_codeml_tasks(
    cds: tuple[SeqRecord, SeqRecord],
    pep: tuple[SeqRecord, SeqRecord],
    align_method: Literal["muscle", "mafft"],
    config: ConfigCODEML | CODEML,
) -> dict:
    config = CODEML.escape(config) if isinstance(config, ConfigCODEML) else config
    with tempfile.TemporaryDirectory(dir="/dev/shm") as tmpdir:
        # Step0: Write records as file.
        SeqIO.write(cds, f"{tmpdir}/cds.fna", "fasta")

        # Step1: pep Alignment.
        SeqIO.write(pep, f"{tmpdir}/pep.faa", "fasta")
        pep_aligned = _alignp(align_method, tmpdir, "pep.faa")

        # Step2: Codon Alignment.
        codon = _codon_align(cds, pep_aligned)

        # Step3: do codeml.
        result = config.preset_dnds(tmpdir, codon)
        result["seq1"] = cds[0].id
        result["seq2"] = cds[1].id
        return result


def calc_codeml(
    cds1: str | list[SeqRecord],
    cds2: str | list[SeqRecord],
    pep1: str | list[SeqRecord],
    pep2: str | list[SeqRecord],
    align_method: Literal["muscle", "mafft"],
    config: ConfigCODEML | CODEML,
    processes_num: int = 2,
) -> list:
    """Calculate the dN and dS values between homologous gene pairs by codeml."""
    codeml_start = time.time()
    # Prepare data.
    total_records = []
    # Load if input is string (file path)
    for records in [cds1, cds2, pep1, pep2]:
        total_records.append(
            [record for record in SeqIO.parse(records, "fasta")]
            if isinstance(records, str)
            else records
        )

    # Packup cds, pep
    cds, pep = list(zip(total_records[0], total_records[1])), list(
        zip(total_records[2], total_records[3])
    )

    # Packup all parameters
    args = [(cds[i], pep[i], align_method, config) for i in range(len(cds))]
    chunks_size = max(1, ceil(len(args) / processes_num))

    with Progress() as progress:
        task = progress.add_task("CODEML running.", total=len(args))

        with Pool(processes=processes_num) as pool:
            async_results = [
                pool.apply_async(
                    _calc_codeml_tasks,
                    a,
                    callback=lambda _: progress.update(
                        task, advance=1, chunks=chunks_size
                    ),
                )
                for a in args
            ]
            results = [ar.get() for ar in async_results]

    used_seconds = time.time() - codeml_start
    console.log(
        f"CODEML pairwise dn/ds calculation finished in {used_seconds//60:.0f} min {used_seconds%60:.1f} sec."
    )
    return results
