#
#
#
#
#


import re

from Bio import SeqIO


def re_chr(
    gff3_file_path: str,
    fna_file_path: str,
    map_file_path: str,
    new_gff3_name: str = "",
    new_fna_name: str = "",
    remove_unmapped: bool = False,
) -> None:
    """
    Replace the chromosome name by map file (in `.csv`), and if needed,
    you can remove other chromosome entries which aren't appear on map.

    `gff3_file_path` `fna_file_path` `map_file_path`  \n
    `new_gff3_name` `new_fna_name`    \n
    `remove_unmapped`   \n
    """

    def switch_dict(map_path: str) -> dict:
        dictionary = {}
        with open(map_path, "r") as map:
            for line in map:
                line_ls = line.strip().split(",")
                dictionary[line_ls[0]] = line_ls[1]
        return dictionary

    d = switch_dict(map_file_path)
    fna_record = []
    gff_record = ""

    # 更改 FNA 里的条目
    for rec in SeqIO.parse(fna_file_path, "fasta"):
        # 改名；若无对应则保持原状
        old_id = rec.id
        try:
            rec.id = d[rec.id]
        except KeyError:
            rec.id = rec.id
        # 添加记录
        if (not remove_unmapped) or (remove_unmapped and (old_id != rec.id)):
            fna_record.append(rec)

    SeqIO.write(fna_record, new_fna_name, "fasta")

    # 更改 GFF 里的条目
    gff_buffer = []
    pattern = re.compile("|".join(re.escape(k) for k in d.keys()))
    with open(gff3_file_path) as gff_original, open(new_gff3_name, "w") as gff_changed:
        for feature_old in gff_original:
            feature_new = pattern.sub(lambda m: d[m.group(0)], feature_old)
            if (not remove_unmapped) or (
                remove_unmapped and (feature_old != feature_old)
            ):
                gff_buffer.append(feature_new)
        gff_changed.write("".join(gff_record))

    print("DONE!")
