import os
import re
import shutil
import string
import subprocess
from functools import partial
from itertools import product
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Any, List, Literal, Optional

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)

tasks_count = 0
completed = 0


def run_single_command(cmd, capture_output=False):
    logger.info(f"Running {' '.join(cmd)}")
    try:
        if capture_output:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return True, result.stdout, result.stderr
        else:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
            return True, None, None
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if hasattr(e, "stderr") else str(e)
        logger.error(
            f"Program {cmd[0] if cmd else 'unknown'} exited abnormally, "
            f"return code {e.returncode}. Error message: {error_msg}"
        )
        return False, None, error_msg


def inc_progress(*args):
    global completed
    completed += 1

    if completed % 10 == 0:
        logger.info(
            f"Progress: {completed} / {tasks_count} ({completed / tasks_count * 100:.2f}%)"
        )


class CMDBuilder:
    def __init__(self, program_name: str, program_path: Optional[Path] = None):
        if (program_path is None) or (
            isinstance(program_path, Path) and not program_path.exists()
        ):
            logger.warning(
                f"Not found or provide the path of executable program {program_name}. Trying search from environment variables."
            )
            # try find it from env-var
            cmdpgpath = shutil.which(program_name)
            if cmdpgpath:
                logger.info(
                    f"Get the path of program {program_name} from environment variables: {cmdpgpath}"
                )
            else:
                logger.exception(f"Program {program_name} not found")
                raise FileNotFoundError
        else:
            logger.info(
                f"Using the provided path of program {program_name}: {program_path}"
            )
            cmdpgpath = program_path

        self.program_name = program_name
        self.program_path = cmdpgpath
        self._cmd = [str(cmdpgpath)]
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

    def add_custom_params(self, *params: Any):
        """
        Add other param(s) if you want.
        """
        self._cmd.extend([str(param) for param in params])

    def run(self, *args, dry_run: bool = False, capture_output: bool = False):
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
                f"Dry-run. This command will be executed: '{' '.join(self.cmd)}'"
            )


class BatchCMDBuilder(CMDBuilder):
    def __init__(
        self,
        program_name: str,
        program_path: Optional[Path] = None,
        substitute_method: Literal["paired", "cross"] = "paired",
        force_to_keep_order: bool = False,
    ):
        super().__init__(program_name, program_path)
        self.substitute_method = substitute_method
        self.substitute_place: List[int] = []
        self.substitute_list: List[List[str]] = []
        self._template_place: List[int] = []
        self._template_strs: List[str] = []
        self.force_to_keep_order = force_to_keep_order

    def add_substitute_param(
        self, param_name: Optional[str], substitute_list: List[Any]
    ):
        if param_name is None:
            self._cmd.append(f"subs {len(self.substitute_place)}")
        else:
            self._cmd.extend([param_name, f"subs {len(self.substitute_place)}"])
        self.substitute_place.append(len(self._cmd) - 1)
        self.substitute_list.append(list(substitute_list))

    def add_substitute_template(self, param_name: Optional[str], template: str):
        """
        Add a parameter whose value is a format-string template referencing
        substitution slots by index.

        In cross / paired mode the template is rendered with the current
        combination's values.  Example::

            cmd.add_substitute_param("--in",  ["X", "Y"])     # slot 0
            cmd.add_substitute_param(None,    ["1", "2"])      # slot 1
            cmd.add_substitute_template("--out", "{0}_{1}.out")

        generates ``--out`` values like ``X_1.out``, ``Y_2.out``, …

        Template fields MUST use explicit integer indices (``{0}``, ``{1}``,
        …) that refer to already-registered substitution slots.  Implicit
        ``{}`` and named fields are rejected.

        **Path-friendly attribute access.**  When a substitution value is a
        :class:`pathlib.Path` (or any object with attributes), you can use
        Python's native format-field syntax to extract components::

            >>> from pathlib import Path
            >>> cmd = BatchCMDBuilder(program_name="example")
            >>> in_files = [Path(\"data/GCF_001.dmnd\"), Path(\"data/GCF_002.dmnd\")]
            >>> cmd.add_substitute_param(\"--in\", in_files)        # slot 0
            >>> cmd.add_substitute_template(\"--out\", \"{0.stem}.out\")  # just the stem

        Common accessors for :class:`~pathlib.PurePath`:
        ``{N.stem}`` (name without suffix), ``{N.name}`` (final component),
        ``{N.parent.name}`` (enclosing directory), ``{N.suffix}``.
        This works for *any* object attribute, not just paths.
        """
        # --- validate referenced slots ------------------------------------
        refs: List[int] = []
        for _literal, field_name, _fmt_spec, _conv in string.Formatter().parse(
            template
        ):
            if field_name is None:
                continue
            if field_name == "":
                raise ValueError(
                    f"Template {template!r} uses implicit '{{}}'. "
                    "Use explicit indices like {0}, {1} to reference "
                    "substitution slots."
                )
            # Accept "0", "0.stem", "0.parent.name", "0[0]", …
            _idx = re.match(r"(\d+)", field_name)
            if _idx is None:
                raise ValueError(
                    f"Template field {field_name!r} in {template!r} does not "
                    "start with a substitution-slot index. "
                    "Use {0}, {1}, … (optionally with .attr or [key] access)."
                ) from None
            refs.append(int(_idx.group(1)))

        if not refs:
            raise ValueError(
                f"Template {template!r} does not reference any substitution "
                "slot. Use add_param() or add_substitute_param() for "
                "non-template parameters."
            )

        max_ref = max(refs)
        if max_ref >= len(self.substitute_list):
            raise IndexError(
                f"Template {template!r} references slot {max_ref}, but only "
                f"{len(self.substitute_list)} substitution(s) have been "
                f"registered (valid: 0–{max(0, len(self.substitute_list) - 1)}). "
                "Call add_substitute_param() first."
            )

        # --- store ---------------------------------------------------------
        marker = f"tmpl {len(self._template_place)}"
        if param_name is None:
            self._cmd.append(marker)
        else:
            self._cmd.extend([param_name, marker])
        self._template_place.append(len(self._cmd) - 1)
        self._template_strs.append(template)

    def _format_templates(self, cmd: List[str], subs_values: List[str]) -> None:
        """Fill every template placeholder in *cmd* with the current
        substitution values (in-place)."""
        for pos, tmpl in zip(self._template_place, self._template_strs):
            cmd[pos] = tmpl.format(*subs_values)

    @property
    def _cmds(self):
        # No substitutions at all — just the base command
        if len(self.substitute_place) == 0:
            return [self._cmd.copy()]

        # Single or multiple substitutions: build the value-set iterator
        if self.substitute_method == "paired":
            # zip stops at the shortest list; single-list case naturally
            # produces the same 1-tuples as the old len==1 branch
            value_sets = zip(*self.substitute_list)
        else:  # cross
            value_sets = product(*self.substitute_list)

        cmds = []
        for values in value_sets:
            new_cmd = self._cmd.copy()
            # fill independent substitution slots
            for pos, val in zip(self.substitute_place, values):
                new_cmd[pos] = str(val)
            # render template slots that reference multiple substitutions
            self._format_templates(new_cmd, list(values))
            cmds.append(new_cmd)

        return cmds

    def __len__(self):
        return len(self.cmd)

    def run(
        self,
        *args,
        processes: int = 8,
        dry_run: bool = False,
        capture_output: bool = False,
    ):
        global tasks_count
        tasks_count = len(self._cmds)
        logger.info(f"Main process ID: {os.getpid()}")
        if not dry_run:
            if processes == 1:
                for cmd in self._cmds:
                    run_single_command(cmd, capture_output)
                    inc_progress()
            else:
                with Pool(processes=processes) as pool:
                    func = partial(run_single_command, capture_output=capture_output)
                    async_results = [
                        pool.apply_async(func, kwds={"cmd": cmd}, callback=inc_progress)
                        for cmd in self._cmds
                    ]
                    results = [ar.get() for ar in async_results]
                    pool.close()
                    pool.join()
                    logger.info("Process completed.")
        else:
            logger.info("Dry-run. These commands will be executed:")
            for cmd in self._cmds:
                logger.info(f"{' '.join(cmd)}")
                inc_progress()
