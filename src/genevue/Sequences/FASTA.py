from pathlib import Path
from typing import List, Optional
from collections import defaultdict

from genevue import console
from genevue import setup_rich_logger
from genevue.External.HMMER3.HMMSearch import HMMSearch
from genevue.utils import deduplicate

logger = setup_rich_logger(__name__, console)


class FASTA:
    def __init__(
        self,
        sequence_path: Path,
    ):
        self.sequence_path = sequence_path
        self.sequence_idx_path = (
            self.sequence_path.parent / f"{self.sequence_path.stem}.idx"
        )

        self.target_ls = []
        self.target_count = len(self.target_ls)
        self.sequence_count = 0
        self.res_sequence_count = 0

    def makeidx(self, idx_path: Optional[Path] = None):
        counter = 0
        sub_counter = 1

        if idx_path is None:
            self.sequence_idx_path = idx_path

        with (
            open(self.sequence_path, "r") as fin,
            open(self.sequence_idx_path, "w") as fout,
        ):
            for line_idx, line in enumerate(fin):
                if line.startswith(">"):
                    sequence_label = line.strip().split()[0].replace(">", "")
                    sequence_line_idx = line_idx
                    if sub_counter == 1:
                        fout.write(f"{sequence_label} {sequence_line_idx}")
                    else:
                        fout.write(
                            f" {sub_counter}\n{sequence_label} {sequence_line_idx}"
                        )
                        sub_counter = 1
                    counter += 1
                else:
                    sub_counter += 1
            fout.write(f" {sub_counter}\n")
        return counter

    def filter(self, target_path: Path, target_type: str, out_path: Path):
        if target_type not in self.support_target_format():
            logger.critical(
                f"Not supported format! Valid format: {", ".join(self.support_target_format())}"
            )
            raise ValueError

        self._router(target_type, target_path)

        logger.info(f"{target_type} sequence target found.")

        if not self.target_ls:
            logger.warning("No target sequences. Abort.")
            return

        if not self.sequence_idx_path.exists():
            logger.info(f"Generating index file for {self.sequence_path}.")
            self.makeidx()

        seq_start_end = {}

        with open(Path(self.sequence_idx_path)) as f:
            for line in f:
                seqid, seqlineidx, seqlinelength = line.strip().split()
                seq_start_end[seqid] = [
                    int(seqlineidx),
                    int(seqlineidx) + int(seqlinelength) - 1,
                ]

        # reordering
        new_target_ls = []
        for seqid in seq_start_end.keys():
            if seqid in self.target_ls:
                new_target_ls.append(seqid)

        if len(new_target_ls) < len(self.target_ls):
            logger.warning(
                f"Some sequence id not in sequence file: {', '.join(set(self.target_ls) - set(new_target_ls))}"
            )

        self.target_ls = new_target_ls

        start_end = [seq_start_end[seqid] for seqid in self.target_ls]

        with open(self.sequence_path) as fin, open(out_path, "w") as fout:
            counter = 0
            for idx, line in enumerate(fin):
                if counter >= len(start_end):
                    break
                start, end = start_end[counter]
                if start <= idx <= end:
                    fout.write(line)
                if idx == end:
                    counter += 1

    @staticmethod
    def support_target_format():
        return ["plain", "blast6", "hmmsearch_noargs"]

    def _router(self, target_type: str, target_path: Path):
        match target_type:
            case "plain":
                self._parse_plain(target_path)
            case "blast6":
                self._parse_blast6(target_path)
            case "hmmsearch_noargs":
                self._parse_hmmsearch_noargs(target_path)
            case _:
                logger.error(f"Unimplemented branch: {target_type}")

    def _parse_plain(self, target_path: Path):
        with open(target_path, "r", encoding="utf-8") as f:
            for target in f:
                target = target.strip()
                self.target_ls.append(target)

        self.target_ls = deduplicate(self.target_ls)

    def _parse_blast6(self, target_path: Path):
        pass

    def _parse_hmmsearch_noargs(self, target_path: Path):
        hmmsearch_inst = HMMSearch(res_path=target_path)
        hmmsearch_inst.parse_result()

        self.target_ls = hmmsearch_inst.result_entries_filter(max_evalue=0.01)

    @property
    def seq_ids(self) -> List[str]:
        res = []
        with open(self.sequence_path) as f:
            for line in f:
                if line.startswith(">"):
                    res.append(line.removeprefix(">").strip().split()[0])
        return res

    @property
    def seq_descriptions(self):
        res = []
        with open(self.sequence_path) as f:
            for line in f:
                if line.startswith(">"):
                    res.append(line.removeprefix(">").strip())
        return res

    @property
    def seq_lengths(self) -> dict[str, int]:
        res = defaultdict(int)
        seq_id: str = ""
        with open(self.sequence_path) as f:
            for line in f:
                if line.startswith(">"):
                    seq_id = line.removeprefix(">").strip().split()[0]
                else:
                    res[seq_id] += len(line.strip())
        return dict(res)

    def export_seq_ids(self, out_path: Path):
        with open(out_path, "w") as f:
            f.write("\n".join(self.seq_ids))

    def seq_rename_by_order(
        self,
        out_path: Path,
        prefix: str = "",
        suffix: str = "",
        same_width: bool = False,
    ):
        if same_width:
            counter = 0
            with open(self.sequence_path) as f:
                for line in f:
                    if line.startswith(">"):
                        counter += 1
            int_width = len(str(counter))
        else:
            int_width = 0

        with open(self.sequence_path) as fin, open(out_path, "w") as fout:
            counter = 1
            for line in fin:
                if line.startswith(">"):
                    line = f">{prefix}{counter:0{int_width}}{suffix} {line.removeprefix(">")}"
                    counter += 1
                fout.write(line)

    def seq_rename_by_table(self, table_path: Path, out_path: Path):
        # fetch
        substitute_d = {}
        with open(table_path) as f:
            for line in f:
                linels = line.strip().split()
                if len(linels) > 2:
                    logger.warning(
                        "Number of columns is greater than 2, the excess will be ignored."
                    )
                    logger.warning(f"Warning happened on {line.strip()}")
                if len(linels) < 2:
                    logger.warning(
                        "Number of columns is lesser than 2, the excess will be ignored."
                    )
                    logger.warning(f"Warning happened on {line.strip()}")
                    continue
                substitute_d[linels[0]] = linels[1]

        with open(self.sequence_path) as fin, open(out_path, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    seq_id = line.removeprefix(">").strip().split()[0]
                    fout.write(
                        ">"
                        " ".join([substitute_d.get(seq_id, ""), line.removeprefix(">")])
                    )
                else:
                    fout.write(line)
