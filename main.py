import pandas as pd

from genevue.Evo.Coll.powercoll import Collinearity
from genevue.Formats.GFF3 import GFF3
from genevue.Plot.Collinearity import MacroCollPlotMultiple
from genevue.Utils.Parse import blast6reader

id1 = "GCA_019202795.1"
id2 = "GCA_019202805.1"

path1 = f"/run/media/hyli360/9cedf66a-dc7d-4c27-bac5-f9f4b7831b33/Genome/workstation/bgff3/{id1}.parquet"
path2 = f"/run/media/hyli360/9cedf66a-dc7d-4c27-bac5-f9f4b7831b33/Genome/workstation/bgff3/{id2}.parquet"
blast = f"/run/media/hyli360/9cedf66a-dc7d-4c27-bac5-f9f4b7831b33/Genome/workstation/BLASTP_RES/{id1}_{id2}.out"

chr1 = [
    "CM032913.1",
    "CM032914.1",
    "CM032915.1",
    "CM032916.1",
    "CM032917.1",
    "CM032918.1",
    "CM032919.1",
    "CM032920.1",
    "CM032921.1",
    "CM032922.1",
    "CM032923.1",
    "CM032924.1",
    "CM032925.1",
]
chr2 = [
    "CM032900.1",
    "CM032901.1",
    "CM032902.1",
    "CM032903.1",
    "CM032904.1",
    "CM032905.1",
    "CM032906.1",
    "CM032907.1",
    "CM032908.1",
    "CM032909.1",
    "CM032910.1",
    "CM032911.1",
    "CM032912.1",
]

gff3a = GFF3()
gff3a.load_from_file(path1)

gff3b = GFF3()
gff3b.load_from_file(path2)

anchor = pd.read_csv(f"{id1}_{id2}.anchor", sep="\t", comment="#")

mcpm = MacroCollPlotMultiple()

mcpm.add_genome(id1, gff3a.db.to_pandas(), 0.2, "green", chrom_filter=chr1)
mcpm.add_genome(id2, gff3b.db.to_pandas(), 0.8, "blue", chrom_filter=chr2)
mcpm.add_collbox(id1, id2, anchor, "grey")
mcpm.add_line(id1, "CM032913.1", 24456, id2, "CM032900.1", 7100479, "red")
mcpm.draw()
mcpm.save("res.png")
