from genevue.GXF.GFF3tools import BlazingGFF3
from pathlib import Path

bgff3 = BlazingGFF3(Path("Pfi_Features.bgff3"))
bgff3.brief_report()
