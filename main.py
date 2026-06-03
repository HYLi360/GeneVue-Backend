from pathlib import Path

basepath = Path(
    "/run/media/hyli360/9cedf66a-dc7d-4c27-bac5-f9f4b7831b331/Genome/Datasets/ncbi_dataset/data/GCA_054771515.1/genomic.gff"
)

from genevue.GXF.GFFtools import GenomicFeatureTree

gft = GenomicFeatureTree(basepath)
print(gft.features["gene-Ancab_007154"].attr_dict)
print(gft.parent2featureid["gene-Ancab_007154"])
ls = []
for id in gft.iter_subfeatures("gene-Ancab_007154"):
    ls.append(gft.features[id])
ls.sort()
for i in ls:
    print(i)
