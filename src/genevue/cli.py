from typing import Callable

import typer

from genevue.Remote import app_remote
from genevue.QC import app_qc
from genevue.sequences import app_sequence
from genevue.Tools.completion import install_completion, uninstall_completion
from genevue.configure import Configure


def _build_config_commands(_: Callable, configure: Configure):
    """
    A sub builder to specially defile `config` command.
    """
    import typer

    app_config = typer.Typer(help=_("View or edit the configuration."))

    @app_config.command(help=_("Show current configuration."))
    def show():
        from rich import print as rprint

        rprint(configure.config)

    @app_config.command(help=_("Reset configuration to defaults."))
    def reset():
        configure.reset()
        print(_("Configuration reset to defaults."))

    @app_config.command(help=_("Set API key for a service provider."))
    def set_apikey(
        provider: str = typer.Option(..., help=_("Service provider (e.g. NCBI)")),
        key: str = typer.Option(..., help=_("API key")),
    ):
        configure.set_apikey(provider, key)
        configure.save()

    @app_config.command(help=_("Set contact email."))
    def set_email(
        email: str = typer.Option(..., help=_("Email address")),
    ):
        configure.set_email(email)

    return app_config


def _build_app(_: Callable, configure: Configure) -> typer.Typer:
    """
    A main app builder.
    """
    # main app
    app = typer.Typer(
        help=_(
            "A easy-to-use, out-of-box toolkit for bioinformatics analysis and processing. "
            "It can run across different platforms seamlessly, and features an elegant, streamlined command-line interface."
        ),
        context_settings={"help_option_names": ["-h", "--help"]},
        no_args_is_help=True,
        add_completion=True,
    )

    # attach sub-module
    app_tools = typer.Typer(help=_("Utility tools"))
    app.add_typer(app_tools, name="tools")

    app.add_typer(app_qc, name="qc")

    pipeline = typer.Typer(help=_(" "))
    app.add_typer(pipeline, name="pipeline")

    app.add_typer(app_sequence, name="sequence")

    app.add_typer(app_remote, name="remote")

    app_config = _build_config_commands(_, configure)
    app.add_typer(app_config, name="config")

    @app.command(
        name="install_completion",
        help=_("Install auto-completion module."),
    )
    def install():
        install_completion()

    @app.command(
        name="uninstall_completion",
        help=_("Uninstall auto-completion module."),
    )
    def uninstall():
        uninstall_completion()

    return app
