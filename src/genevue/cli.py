from typing import Literal, Optional
from pathlib import Path

import rich.json
import typer

from genevue import console
from genevue import local_configure
from genevue.Tools.completion import install_completion, uninstall_completion
from genevue.l10n import _
from genevue.Pipelines.GeneFamilySearch import Pipeline_GeneFamilySearch

from genevue.QC import app_qc

app = typer.Typer(
    help=_(
        "A easy-to-use, out-of-box toolkit for bioinformatics analysis and processing. "
        "It can run across different platforms seamlessly, and features an elegant, streamlined command-line interface."
    ),
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
    add_completion=False,
)

app_tools = typer.Typer(help=_("这里是一些文字"))
app.add_typer(app_tools, name="tools")


@app_tools.command()
def fasta_renamer(
    mode: Literal["make", "rename"] = typer.Option(None, "--mode", "-m", help=_(" ")),
    input_path: str = typer.Option(None, "--in", "-i", help=_(" ")),
):
    pass


# ---QC----------------------------------------------------------------------------------------------
app.add_typer(app_qc, name="qc")

# ---pipelines----------------------------------------------------------------------------------------
pipeline = typer.Typer(help=_(" "))
app.add_typer(pipeline, name="pipeline")


@pipeline.command(
    help=_("Search gene family by using BLASTp and HMMSEARCH."),
)
def search_gene_family(
    pep_path: str = typer.Option(None, "--pep", "-p", help=_("PEP sequences path.")),
    blastp_probe_path: Optional[str] = typer.Option(
        None, "--blastp-probe", help=_(" ")
    ),
    blastp_pdb_path: str = typer.Option(None, "--blastp-pdb"),
    blastp_res_path: str = typer.Option(None, "--blastp-res"),
    blastp_max_evalue: float = typer.Option(0.1, "--blastp-max-evalue"),
    blastp_min_bitscore: float = typer.Option(0, "--blastp-min-bitscore"),
    phmm_probe_path: str = typer.Option(None, "--phmm-probe"),
    hmm_max_evalue: float = typer.Option(0.1, "--hmm-max-evalue"),
    hmm_res_path: str = typer.Option(None, "--hmm-res"),
    record_output_path: str = typer.Option(None, "--out", "-o"),
):
    if pep_path is None:
        pep_path = typer.prompt("Set the PEP sequences path. ", prompt_suffix="@ ")

    if blastp_probe_path is None:
        blastp_probe_path = typer.prompt(
            "Set the BLASTp probe path. ", prompt_suffix="@ "
        )
    if phmm_probe_path is None:
        phmm_probe_path = typer.prompt(
            "Set the HMMSEARCH probe path. ", prompt_suffix="@ "
        )
    if blastp_res_path is None:
        blastp_res_path = typer.prompt(
            "blastp res path? ",
            prompt_suffix="@ ",
            default="blastp.res",
            show_default=True,
        )

    if hmm_res_path is None:
        hmm_res_path = typer.prompt(
            "hmm res path? ",
            prompt_suffix="@ ",
            default="hmmsearch.res",
            show_default=True,
        )

    if record_output_path is None:
        record_output_path = typer.prompt("record output? ", prompt_suffix="@ ")

    if_change_max_evalue = typer.confirm(
        "Change the max-evalue or min-bitscore? ",
        prompt_suffix="> ",
        default=False,
        show_default=True,
    )
    if if_change_max_evalue:
        blastp_max_evalue = typer.prompt(
            "Max BLASTp evalue? ", prompt_suffix="> ", default=0.1, show_default=True
        )
        blastp_min_bitscore = typer.prompt(
            "Min BLASTp bitscore? ", prompt_suffix="> ", default=0, show_default=True
        )
        hmm_max_evalue = typer.prompt(
            "Max HMMSEARCH evalue? ", prompt_suffix="> ", default=0.1, show_default=True
        )

    Pipeline_GeneFamilySearch.execute(
        input_data={
            "pep_path": pep_path,
            "blastp_probe_path": blastp_probe_path,
            "blastp_pdb_path": blastp_pdb_path,
            "blastp_res_path": blastp_res_path,
            "blastp_max_evalue": blastp_max_evalue,
            "blastp_min_bitscore": blastp_min_bitscore,
            "phmm_probe_path": phmm_probe_path,
            "hmm_max_evalue": hmm_max_evalue,
            "hmm_res_path": hmm_res_path,
            "record_output_path": record_output_path,
        }
    )


# ---Remote-------------------------------------------------------------------------------------------
app_remote = typer.Typer(help=_("Some useful remote or download tools."))
app.add_typer(app_remote, name="remote")


# ---Config-------------------------------------------------------------------------------------------
app_config = typer.Typer(help=_("View or edit the configuration."))
app.add_typer(app_config, name="config")


@app_config.command()
def path():
    pass


@app_config.command()
def reset():
    is_confirmed = typer.confirm(
        "Do you want to reset the configuration? ",
        prompt_suffix="> ",
    )
    if is_confirmed:
        local_configure.reset()
        console.print(f"Configuration reset successfully!")


@app_config.command()
def disp():
    console.print(local_configure.config)


@app_config.command()
def email(new_email: str):
    if not new_email:
        console.warn("You didn't input anything. Abort.")
    else:
        console.info(f"Old e-mail address: {local_configure.get_email()}")
        local_configure.set_email(new_email)
        console.info(f"New e-mail address: {new_email}")


# ---Install (or uninstall) auto-completion module----------------------------------------------------
@app.command(
    help=_("Install auto-completion module."),
)
def install():
    install_completion()


@app.command(
    help=_("Uninstall auto-completion module."),
)
def uninstall():
    uninstall_completion()
