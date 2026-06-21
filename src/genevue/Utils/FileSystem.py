import gzip
from pathlib import Path
from types import NoneType
from typing import Optional, List

import bgzip

from genevue import setup_rich_logger, console

logger = setup_rich_logger(__name__, console)

MAGIC_NUMBER_HEAD = {
    "gz": (3, b"\x1f\x8b\x08"),
    "gzip": (3, b"\x1f\x8b\x08"),
    "gff3": (15, b"##gff-version 3"),
    "parquet": (4, b"\x50\x41\x52\x31"),
}


def file_renamer(filename_list: list[str], prefix: str = "", suffix: str = ""):
    return ["".join([prefix, filename, suffix]) for filename in filename_list]


def echo_my_cwd():
    console.print(f"Your current directory: {Path.cwd()}")


def iter_path(
    target_path: Path = Path.cwd(),
    output: Optional[Path] = None,
    include_folder: bool = False,
    norec: bool = True,
) -> List[Path]:
    if not target_path:
        target_path = Path.cwd()

    if not target_path.resolve().exists():
        console.print(f"Path {target_path.resolve()} not exists. [red]abort.[/red]")
        return []

    file_handler = open(output, "w") if isinstance(output, Path) else None

    if not isinstance(file_handler, NoneType):
        if norec:
            for subpath in target_path.iterdir():
                if subpath.is_dir() and not include_folder:
                    continue
                file_handler.write(f"{subpath}\n")
        else:
            if include_folder:
                for root, dirs, files in target_path.walk():
                    for name in dirs:
                        file_handler.write(f"{root / name}\n")
                    for name in files:
                        file_handler.write(f"{root / name}\n")
            else:
                for root, dirs, files in target_path.walk():
                    for name in files:
                        file_handler.write(f"{root / name}\n")
        return []
    else:
        res: List[Path] = []
        if norec:
            for subpath in target_path.iterdir():
                if subpath.is_dir() and not include_folder:
                    continue
                res.append(subpath)
        else:
            if include_folder:
                for root, dirs, files in target_path.walk():
                    for name in dirs:
                        res.append(root / name)
                    for name in files:
                        res.append(root / name)
            else:
                for root, dirs, files in target_path.walk():
                    for name in files:
                        res.append(root / name)
        res.sort()
        return res


def default_filename(
    file_name_ls: List[Path],
    prefix: str = "",
    suffix: str = "",
    exname: str = "",
) -> List[Path]:
    new_paths = []
    for path in file_name_ls:
        new_paths.append(
            path.parent
            / f"{prefix}{path.stem}{suffix}{f".{exname}" if exname and path.suffix else path.suffix}"
        )
    return new_paths


def pair_filename(
    input_dir: Path,
    output_dir: Path,
    prefix: str = "",
    suffix: str = "",
    exname: str = "",
):
    in_file_name_stem_ls = [file.stem for file in iter_path(input_dir)]
    in_file_name_stem_ls.sort()

    exname = exname.replace(".", "")

    output_file_name_ls = [
        output_dir / f"{prefix}{file_name}{suffix}{f".{exname}" if exname else ""}"
        for file_name in in_file_name_stem_ls
    ]
    output_file_name_ls.sort()

    return iter_path(input_dir), output_file_name_ls


def check_filetype(file_path: Path, target_type: str) -> bool:
    leng, byte = MAGIC_NUMBER_HEAD.get(target_type, (0, b""))
    if leng == 0:
        print(f"Not a valid format: {target_type}")
        return False
    with open(file_path, "rb") as f:
        if byte == f.read(leng):
            return True
    return False


def gzip2bgzip(gzip_path: Path, bgzip_path: Path, chunk_size: int = 4096):
    if check_filetype(gzip_path, "gz"):
        with gzip.open(gzip_path, "rb") as gz_f:
            with open(bgzip_path, "wb") as bgz_f:
                bgz_writer = bgzip.BGZipWriter(bgz_f)

                while True:
                    chunk = gz_f.read(chunk_size)
                    if not chunk:
                        break
                    bgz_writer.write(chunk)

                bgz_writer.close()
            console.print("[green]Transform completed!")
    else:
        console.print(f"{gzip_path} is not a valid gzip file! abort.")
