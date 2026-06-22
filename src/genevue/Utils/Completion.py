"""
This code is ported from Typer's `_completion_classes.py` and `_completion_shared.py`,
but has been slightly streamlined to fit into a single script file.

It still supports common Linux shells (Bash, Zsh, and Fish) and Windows PowerShell.

Thank you, <https://tiangolo.com/>!
"""

import os
from pathlib import Path

import shellingham

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)

COMPLETION_SCRIPT_BASH = """
_genevue_completion() {
    local IFS=$'\n'
    COMPREPLY=( $( env COMP_WORDS="${COMP_WORDS[*]}" \\
                   COMP_CWORD=$COMP_CWORD \\
                   _GENEVUE_COMPLETE=complete_bash $1 ) )
    return 0
}

complete -o default -F _genevue_completion genevue
"""

COMPLETION_SCRIPT_ZSH = """
#compdef genevue

_genevue_completion() {
  eval $(env _TYPER_COMPLETE_ARGS="${words[1,$CURRENT]}" _GENEVUE_COMPLETE=complete_zsh genevue)
}

compdef _genevue_completion genevue
"""

COMPLETION_SCRIPT_FISH = 'complete --command genevue --no-files --arguments "(env _GENEVUE_COMPLETE=complete_fish _TYPER_COMPLETE_FISH_ACTION=get-args _TYPER_COMPLETE_ARGS=(commandline -cp) genevue)" --condition "env _GENEVUE_COMPLETE=complete_fish _TYPER_COMPLETE_FISH_ACTION=is-args _TYPER_COMPLETE_ARGS=(commandline -cp) genevue"'

COMPLETION_SCRIPT_POWER_SHELL = """
Import-Module PSReadLine
Set-PSReadLineKeyHandler -Chord Tab -Function MenuComplete
$scriptblock = {
    param($wordToComplete, $commandAst, $cursorPosition)
    $Env:_GENEVUE_COMPLETE = "complete_powershell"
    $Env:_TYPER_COMPLETE_ARGS = $commandAst.ToString()
    $Env:_TYPER_COMPLETE_WORD_TO_COMPLETE = $wordToComplete
    genevue | ForEach-Object {
        $commandArray = $_ -Split ":::"
        $command = $commandArray[0]
        $helpString = $commandArray[1]
        [System.Management.Automation.CompletionResult]::new(
            $command, $command, 'ParameterValue', $helpString)
    }
    $Env:_GENEVUE_COMPLETE = ""
    $Env:_TYPER_COMPLETE_ARGS = ""
    $Env:_TYPER_COMPLETE_WORD_TO_COMPLETE = ""
}
Register-ArgumentCompleter -Native -CommandName genevue -ScriptBlock $scriptblock
"""


def _which_shell() -> str:
    # Use `shellingham` to check the type of shell environment
    name: str = shellingham.detect_shell()[0].lower()
    return name


def _get_script(shell_name: str) -> str:
    script_dict = {
        "bash": COMPLETION_SCRIPT_BASH,
        "zsh": COMPLETION_SCRIPT_ZSH,
        "fish": COMPLETION_SCRIPT_FISH,
        "powershell": COMPLETION_SCRIPT_POWER_SHELL,
        "pwsh": COMPLETION_SCRIPT_POWER_SHELL,
    }

    return script_dict.get(shell_name, "")


def install_completion() -> None:
    shell_name = _which_shell()
    script = _get_script(shell_name)
    logger.info(f"Shell type '{shell_name}' detected.")

    def bash():
        completion_script_path = Path.home() / ".bash_completions" / "genevue.sh"
        completion_script_path.parent.mkdir(exist_ok=True, parents=True)
        rc_path = Path.home() / ".bashrc"
        rc_path.parent.mkdir(exist_ok=True, parents=True)

        # Write the script in f".bash_completions/genevue.sh"
        logger.info(f"Write auto-completion module on {completion_script_path}.")
        with open(completion_script_path, "w") as completion_script_handler:
            completion_script_handler.write(script)

        # Append a "source" line after the ".bashrc"
        # Search "source" line before real append
        with open(rc_path, "r") as rc_handler_r:
            for line in rc_handler_r:
                if line == f"source '{completion_script_path}'\n":
                    logger.info(".bashrc has been configured before install.")
                    return

        logger.info(f"Setting user's .bashrc configures. on {rc_path}.")
        with open(rc_path, "a") as rc_handler_w:
            rc_handler_w.write(
                "# ---This script was appended by GeneVue---\n"
                f"source '{completion_script_path}'\n"
                "# ---Run 'genevue uninstall' to purge it and associated script---\n"
            )

    def zsh():
        pass

    def fish():
        pass

    def pwsh():
        pass

    install_instructions_dict = {
        "bash": bash,
    }

    install_instructions_dict.get(shell_name)()
    logger.info(
        "Auto-complete module installed. Changes will take effect when you restart the shell.\n"
        "Run 'genevue uninstall' to remove this module."
    )


def uninstall_completion() -> None:
    shell_name = _which_shell()

    def bash():
        completion_script_path = Path.home() / ".bash_completions" / "genevue.sh"
        rc_path = Path.home() / ".bashrc"

        # Delete ".bash_completions/genevue.sh"
        logger.info(f"Removing {completion_script_path}.")
        try:
            os.remove(completion_script_path)
        except FileNotFoundError:
            console.warn(f"Not found the script {completion_script_path}. Skipping.")

        # Remove the "source" line
        logger.info(f"Remove the source line in {rc_path}.")
        new_bashrc_content = []
        with open(rc_path, "r") as rc_handler_r:
            for line in rc_handler_r:
                if line not in (
                    "# ---This script was appended by GeneVue---\n",
                    f"source '{completion_script_path}'\n",
                    "# ---Run 'genevue uninstall' to purge it and associated script---\n",
                ):
                    new_bashrc_content.append(line)

        with open(rc_path, "w") as rc_handler_w:
            rc_handler_w.write("".join(new_bashrc_content))

    uninstall_instructions_dict = {
        "bash": bash,
    }

    uninstall_instructions_dict.get(shell_name)()
    logger.info("Auto-complete module uninstalled.")
