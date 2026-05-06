import os
import sys
import typer
from typing import Optional, Literal

from genevue import console
from genevue.l10n import _
from genevue.QC.BUSCO import __version__
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
from genevue.QC.BUSCO.BuscoDownloadManager import BuscoDownloadManager
from genevue.QC.BUSCO.BuscoConfig import (
    PseudoConfigForDownload,
    PseudoConfigForPlot,
    BaseConfig,
)
from genevue.QC.BUSCO.ConfigManager import BuscoConfigManager

app_qc = typer.Typer(help=" ")


def _build_params(
    input_sequence_file,
    output_dir,
    process_mode,
    lineage_dataset,
    use_augustus,
    augustus_parameters,
    augustus_species,
    is_auto_lineage,
    is_auto_lineage_euk,
    is_auto_lineage_prok,
    threads_num,
    config_file_path,
    config_break_num,
    dataset_version,
    download,
    download_base_url,
    download_path,
    blast_max_evalue,
    force_rewriting,
    limited,
    is_augustus_long,
    is_using_metaeuk,
    metaeuk_parameters,
    metaeuk_rerun_parameters,
    is_using_miniprot,
    is_skip_bbtools,
    is_offline,
    is_not_telemetry,
    out_path,
    plot,
    is_plot_percentages,
    is_quiet,
    is_restart,
    is_scaffold_composition,
    is_tar_output,
):
    """Map Typer parameter names to the internal key names expected by BuscoConfigManager."""
    params = {}

    if input_sequence_file is not None:
        params["in"] = input_sequence_file
    if output_dir is not None:
        params["out"] = output_dir
    if process_mode is not None:
        params["mode"] = process_mode
    if lineage_dataset is not None:
        params["lineage_dataset"] = lineage_dataset
    if use_augustus:
        params["use_augustus"] = True
    if augustus_parameters is not None:
        params["augustus_parameters"] = augustus_parameters
    if augustus_species is not None:
        params["augustus_species"] = augustus_species
    if is_auto_lineage:
        params["auto-lineage"] = True
    if is_auto_lineage_euk:
        params["auto-lineage-euk"] = True
    if is_auto_lineage_prok:
        params["auto-lineage-prok"] = True
    if threads_num is not None:
        params["cpu"] = threads_num
    if config_file_path is not None:
        params["config_file"] = config_file_path
    if config_break_num is not None:
        params["contig_break"] = config_break_num
    if dataset_version is not None:
        params["datasets_version"] = dataset_version
    if download is not None:
        params["download"] = download
    if download_base_url is not None:
        params["download_base_url"] = download_base_url
    if download_path is not None:
        params["download_path"] = download_path
    if force_rewriting:
        params["force"] = True
    if is_augustus_long:
        params["long"] = True
    if is_using_metaeuk:
        params["use_metaeuk"] = True
    if metaeuk_parameters is not None:
        params["metaeuk_parameters"] = metaeuk_parameters
    if metaeuk_rerun_parameters is not None:
        params["metaeuk_rerun_parameters"] = metaeuk_rerun_parameters
    if is_using_miniprot:
        params["use_miniprot"] = True
    if is_skip_bbtools:
        params["skip_bbtools"] = True
    if is_offline:
        params["offline"] = True
    if is_not_telemetry:
        params["opt-out-run-stats"] = True
    if out_path is not None:
        params["out_path"] = out_path
    if plot is not None:
        params["plot"] = plot
    if is_plot_percentages:
        params["plot_percentages"] = True
    if is_quiet:
        params["quiet"] = True
    if is_restart:
        params["restart"] = True
    if is_scaffold_composition:
        params["scaffold_composition"] = True
    if is_tar_output:
        params["tar"] = True

    if blast_max_evalue is not None and (
        (process_mode != "prot") and (process_mode != "protein")
    ):
        params["evalue"] = blast_max_evalue
    else:
        params["evalue"] = None
    if limited is not None and (
        (process_mode != "prot") and (process_mode != "protein")
    ):
        params["limit"] = limited
    else:
        params["limit"] = None

    return params


def _cleanup_log():
    """Safely close any lingering BUSCO file handler (no-op if none active)."""
    BuscoLogger.close_file_handler()
    try:
        os.remove("busco_{}.log".format(BuscoLogger.pid))
    except FileNotFoundError:
        pass


@app_qc.command(
    help="BUSCO v6.0 for GeneVue.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
    name="busco",
)
def busco(
    input_sequence_file: Optional[str] = typer.Option(
        None,
        "-i",
        "--in",
        metavar="SEQUENCE_FILE",
        help=_(
            "Input sequence file in FASTA format. "
            "Can be an assembled genome or transcriptome (DNA,) or protein sequences from an annotated gene set. "
            "Also possible to use a path to a directory containing multiple input files."
        ),
    ),
    output_dir: Optional[str] = typer.Option(
        None,
        "-o",
        "--out",
        metavar="OUTPUT",
        help=_(
            "Give your analysis run a recognisable short name. "
            "Output folders and files will be labelled with this name. "
            "The path to the output folder is set with --out-path."
        ),
    ),
    process_mode: Optional[
        Literal["geno", "genome", "tran", "transcriptome", "prot", "proteins"]
    ] = typer.Option(
        None,
        "-m",
        "--mode",
        help="Specify which BUSCO analysis mode to run. "
        "There are three valid modes: [green]geno[/green] or [green]genome[/green], for genome assemblies (DNA;) [green]"
        "tran[/green] or [green]transcriptome[/green], for transcriptome assemblies (DNA;) [green]prot[/green] or [green]"
        "proteins[/green], for annotated gene sets (protein.)",
    ),
    lineage_dataset: Optional[str] = typer.Option(
        None,
        "-l",
        "--lineage-dataset",
        metavar="LINEAGE",
        help="Specify the name of the BUSCO lineage to be used.",
    ),
    use_augustus: bool = typer.Option(
        False,
        "--use-augustus",
        is_flag=True,
        help="Use augustus gene predictor for eukaryote runs",
    ),
    augustus_parameters: Optional[str] = typer.Option(
        None,
        "--augustus_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        help="Pass additional arguments to Augustus. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    ),
    augustus_species: Optional[str] = typer.Option(
        None,
        "--augustus_species",
        help="Specify a species for Augustus training.",
    ),
    is_auto_lineage: bool = typer.Option(
        False,
        "--auto-lineage",
        is_flag=True,
        help="Run auto-lineage to find optimum lineage path.",
    ),
    is_auto_lineage_euk: bool = typer.Option(
        False,
        "--auto-lineage-euk",
        is_flag=True,
        help="Run auto-placement just on eukaryote tree to find optimum lineage path",
    ),
    is_auto_lineage_prok: bool = typer.Option(
        False,
        "--auto-lineage-prok",
        is_flag=True,
        help="Run auto-lineage just on non-eukaryote trees to find optimum lineage path",
    ),
    threads_num: Optional[int] = typer.Option(
        None,
        "-c",
        "--cpu",
        metavar="N",
        help="Specify the number of threads/cores to use.",
    ),
    config_file_path: Optional[str] = typer.Option(
        None, "--config", help="Provide a config file."
    ),
    config_break_num: Optional[int] = typer.Option(
        None,
        "--contig_break",
        metavar="n",
        help="Number of contiguous Ns to signify a break between contigs.",
    ),
    dataset_version: Optional[Literal["odb10", "odb12"]] = typer.Option(
        None, "--datasets_version", help="Specify the version of BUSCO datasets"
    ),
    download: Optional[str] = typer.Option(
        None,
        "--download",
        help='Download dataset. Possible values are a specific dataset name, "all", "prokaryota", "eukaryota", '
        'or "virus". If used together with other command line arguments, make sure to place this last.',
    ),
    download_base_url: Optional[str] = typer.Option(
        None,
        "--download_base_url",
        help="Set the url to the remote BUSCO dataset location",
    ),
    download_path: Optional[str] = typer.Option(
        None,
        "--download_path",
        help="Specify local filepath for storing BUSCO dataset downloads",
    ),
    blast_max_evalue: float = typer.Option(
        BaseConfig.BLAST_ARGS["evalue"],
        "-e",
        "--evalue",
        help="E-value cutoff for BLAST searches. Allowed formats, 0.001 or 1e-03",
    ),
    force_rewriting: bool = typer.Option(
        False,
        "-f",
        "--force",
        is_flag=True,
        help="Force rewriting of existing files. Must be used when output files with the provided name already exist.",
    ),
    limited: Optional[int] = typer.Option(
        BaseConfig.BBTOOLS_ARGS["contig_break"],
        "--limited",
        help="How many candidate regions (contig or transcript) to consider per BUSCO",
    ),
    is_list_datasets: bool = typer.Option(
        False,
        "--list-datasets",
        is_flag=True,
        help="Print the list of available BUSCO datasets",
    ),
    is_augustus_long: bool = typer.Option(
        False,
        "--long",
        is_flag=True,
        help="Optimization Augustus self-training mode (Default: Off); adds considerably to the run "
        "time, but can improve results for some non-model organisms",
    ),
    is_using_metaeuk: bool = typer.Option(
        False, "--metaeuk", is_flag=True, help="Use Metaeuk gene predictor"
    ),
    metaeuk_parameters: Optional[str] = typer.Option(
        None,
        "--metaeuk_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        help="Pass additional arguments to Metaeuk for the first run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    ),
    metaeuk_rerun_parameters: Optional[str] = typer.Option(
        None,
        "--metaeuk_rerun_parameters",
        metavar='"--PARAM1=VALUE1,--PARAM2=VALUE2"',
        help="Pass additional arguments to Metaeuk for the second run. All arguments should be contained within a "
        "single string with no white space, with each argument separated by a comma.",
    ),
    is_using_miniprot: bool = typer.Option(
        False,
        "--miniprot",
        is_flag=True,
        help="Use Miniprot gene predictor",
    ),
    is_skip_bbtools: bool = typer.Option(
        False,
        "--skip_bbtools",
        is_flag=True,
        help="Skip BBTools for assembly statistics",
    ),
    is_offline: bool = typer.Option(
        False,
        "--offline",
        is_flag=True,
        help="To indicate that BUSCO cannot attempt to download files",
    ),
    is_not_telemetry: bool = typer.Option(
        False,
        "--opt-out-run-stats",
        is_flag=True,
        help="Opt out of data collection. Information on the data collected is available in the user guide.",
    ),
    out_path: Optional[str] = typer.Option(
        None,
        "--out_path",
        metavar="OUTPUT_PATH",
        help="Optional location for results folder, excluding results folder name. "
        "Default is current working directory.",
    ),
    plot: Optional[str] = typer.Option(
        None,
        "--plot",
        metavar="WORKING_DIRECTORY",
        help="Generate a BUSCO summary plot for all short summary files in the given working directory.",
    ),
    is_plot_percentages: bool = typer.Option(
        False,
        "--plot_percentages",
        is_flag=True,
        help="Plot the percentages of BUSCOs instead of the number of BUSCOs. To be used as an option with --plot.",
    ),
    is_quiet: bool = typer.Option(
        False,
        "-q",
        "--quiet",
        is_flag=True,
        help="Disable the info logs, displays only errors",
    ),
    is_restart: bool = typer.Option(
        False,
        "-r",
        "--restart",
        is_flag=True,
        help="Continue a run that had already partially completed.",
    ),
    is_scaffold_composition: bool = typer.Option(
        False,
        "--scaffold_composition",
        is_flag=True,
        help="Writes ACGTN content per scaffold to a file scaffold_composition.txt",
    ),
    is_tar_output: bool = typer.Option(
        False,
        "--tar",
        is_flag=True,
        help="Compress some subdirectories with many files to save space",
    ),
    is_version: bool = typer.Option(
        False, "-v", "--version", is_flag=True, help="Show this version and exit"
    ),
):
    params = _build_params(
        input_sequence_file,
        output_dir,
        process_mode,
        lineage_dataset,
        use_augustus,
        augustus_parameters,
        augustus_species,
        is_auto_lineage,
        is_auto_lineage_euk,
        is_auto_lineage_prok,
        threads_num,
        config_file_path,
        config_break_num,
        dataset_version,
        download,
        download_base_url,
        download_path,
        blast_max_evalue,
        force_rewriting,
        limited,
        is_augustus_long,
        is_using_metaeuk,
        metaeuk_parameters,
        metaeuk_rerun_parameters,
        is_using_miniprot,
        is_skip_bbtools,
        is_offline,
        is_not_telemetry,
        out_path,
        plot,
        is_plot_percentages,
        is_quiet,
        is_restart,
        is_scaffold_composition,
        is_tar_output,
    )

    if is_version:
        console.print(
            f"\n[blue]BUSCO {__version__}[/blue] for [green]GeneVue[/green].\n"
            "BUSCO: the Benchmarking Universal Single-Copy Ortholog assessment tool.\n"
            "Need help? Check https://busco.ezlab.org/busco_userguide.html\n"
            "Already used BUSCO? Check https://gitlab.com/ezlab/busco#how-to-cite-busco\n"
        )
        sys.exit(0)

    if is_quiet:
        BuscoLogger.quiet = True

    # --list-datasets: replicate ListLineagesAction logic
    if is_list_datasets:
        try:
            config_manager = BuscoConfigManager(params)
            config = PseudoConfigForDownload(
                config_manager.config_file, config_manager.params
            )
            config.load()
            version = dataset_version or "odb12"
            if version == "odb10":
                lineages_list_file = config.downloader.get(
                    "lineages_list_odb10.txt", "information"
                )
            else:
                lineages_list_file = config.downloader.get(
                    "lineages_list.txt", "information"
                )
            with open(lineages_list_file, "r") as f:
                print("".join(f.readlines()))
        except SystemExit as se:
            BuscoLogger.get_logger(__name__).error(se)
            raise SystemExit(1)
        sys.exit(0)

    # --plot: replicate GeneratePlotAction logic
    if plot is not None:
        try:
            from genevue.QC.BUSCO.BuscoPlot import Plot

            config_manager = BuscoConfigManager(params)
            config = PseudoConfigForPlot(
                config_manager.config_file, config_manager.params
            )
            config.load()
            p = Plot(config, plot)
            p.load_data()
            p.generate_plot()
        except SystemExit as se:
            BuscoLogger.get_logger(__name__).error(se)
            raise SystemExit(1)
        sys.exit(0)

    # --download: replicate DirectDownload logic
    if download is not None:
        try:
            config_manager = BuscoConfigManager(params)
            config = PseudoConfigForDownload(
                config_manager.config_file, config_manager.params
            )
            config.load()
            bdm = BuscoDownloadManager(config)
            _download_datasets(bdm, download)
        except SystemExit as se:
            BuscoLogger.get_logger(__name__).error(se)
            raise SystemExit(1)
        sys.exit(0)

    # Main BUSCO run
    from genevue.QC.BUSCO.run_BUSCO import BuscoMaster

    busco_run = BuscoMaster(params)
    busco_run.run()


def _download_datasets(bdm, item):
    """Replicate DirectDownload.get() logic for a single dataset item."""
    if item in ("all", "all_odb12"):
        files_to_get = [i for i in bdm.version_files.index if "odb12" in i]
    elif item in ("prokaryota", "all_prokaryota_odb12"):
        files_to_get = [
            i
            for i in bdm.version_files.index
            if (bdm.version_files.loc[i]["domain"] == "Prokaryota") and ("odb12" in i)
        ]
    elif item in ("eukaryota", "all_eukaryota_odb12"):
        files_to_get = [
            i
            for i in bdm.version_files.index
            if (bdm.version_files.loc[i]["domain"] == "Eukaryota") and ("odb12" in i)
        ]
    elif item == "virus":
        files_to_get = [
            i
            for i in bdm.version_files.index
            if bdm.version_files.loc[i]["domain"] == "Virus"
        ]
    elif item == "all_odb10":
        files_to_get = [i for i in bdm.version_files.index if "odb10" in i]
    elif item == "all_prokaryota_odb10":
        files_to_get = [
            i
            for i in bdm.version_files.index
            if (bdm.version_files.loc[i]["domain"] == "Prokaryota") and ("odb10" in i)
        ]
    elif item == "all_eukaryota_odb10":
        files_to_get = [
            i
            for i in bdm.version_files.index
            if (bdm.version_files.loc[i]["domain"] == "Eukaryota") and ("odb10" in i)
        ]
    elif isinstance(item, str) and item in bdm.version_files.index:
        files_to_get = [item]
    else:
        BuscoLogger.get_logger(__name__).error(
            "{} is not a recognized option".format(item)
        )
        return

    filetypes = [bdm.version_files.loc[f]["type"] for f in files_to_get]
    for f, filename in enumerate(files_to_get):
        bdm.get(filename, filetypes[f])
