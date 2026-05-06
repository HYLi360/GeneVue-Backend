from pathlib import Path
import os, gzip, pysam
from genevue.utils import header_detector
from genevue import FileNotExistsError, console
from genevue.l10n import _


def bgzip_recompressor(ingz: str | Path, outbgz: str | Path):
    ingz = Path(ingz).resolve()
    outbgz = Path(outbgz).resolve()
    if not ingz.exists():
        raise FileNotExistsError(ingz)

    # check the header of input
    if header_detector(ingz, "bgz"):
        console.log(_("There is no need to recompress!"))
        console.log(_("Abort."))
        return 0

    with gzip.open(ingz, "rb") as fin, pysam.BGZFile(str(outbgz), "wb", None) as fout:
        for chunk in iter(lambda: fin.read(1024 * 1024), b""):
            fout.write(chunk)

    console.log(_("gz->bgz done!"))


def filepath_collector_bettor(
    target_path: str | Path = "",
    path_scheme: str = "",
):
    target_path = (
        Path(target_path).resolve() if isinstance(target_path, str) else target_path
    )
    try:
        open(target_path)
    except:
        console.exception(FileNotExistsError(target_path))
