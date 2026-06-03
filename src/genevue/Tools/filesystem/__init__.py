from types import NoneType

from genevue import setup_rich_logger, console
from pathlib import Path

from typing import Optional, List

logger = setup_rich_logger(__name__, console)


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
