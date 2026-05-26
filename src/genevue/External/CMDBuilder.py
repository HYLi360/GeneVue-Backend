from pathlib import Path
from typing import Any

import subprocess
import shutil

from genevue import setup_rich_logger, console

logger = setup_rich_logger(__name__, console)


class CMDBuilder:
    def __init__(self, program_name: str, program_path: Path):
        self.program_name = program_name
        self.program_path = program_path
        self._cmd = [str(program_path)]
        self.result = None

    @property
    def cmd(self):
        return self._cmd

    def add_param(self, param_name: str, param_value: Any):
        """
        Add a 'key-value' type param for this command. `-` or `--` may needed.
        """
        self._cmd.extend([param_name, str(param_value)])

    def add_flag(self, flag_name: str, flag_bool: bool = True):
        """
        Add a 'key' type param for this command. `-` or `--` may needed.
        """
        if flag_bool:
            self._cmd.append(flag_name)

    def run(self, *args, dry_run: bool = False, capture_output: bool = False):
        if not self.program_path.exists():
            # try find it from env-var
            cmdpgpath = shutil.which(self.program_name)
            if cmdpgpath:
                logger.info(
                    f"Get the program path from environment varient: {cmdpgpath}"
                )
            else:
                logger.exception(f"Program not found on {self.program_path}")
                raise FileNotFoundError
        else:
            logger.info(
                f"Get the program path by configure file of genevue: {self.program_path}"
            )

        if args:
            self._cmd.extend([str(arg) for arg in args])

        if not dry_run:
            try:
                logger.info(f"Run the external command: {' '.join(self.cmd)}")
                self.result = subprocess.run(self.cmd, capture_output)
            except subprocess.CalledProcessError as e:
                logger.exception(
                    f"External program {self.program_name} raised a exception with return code {e.returncode}\n"
                    "Std-Err output:\n"
                    f"{e.stderr}"
                )
        else:
            logger.info(
                f"Dry-run. This command will be executed: `{' '.join(self.cmd)}`"
            )
