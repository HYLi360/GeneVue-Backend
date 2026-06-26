#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

from typing import Annotated, Callable, Literal

import typer

from genevue import __full_version__, console
from genevue.configure import Configure
from genevue.Diagnosis import Diagnosis
from genevue.Preprocessing import app_preprocessing
from genevue.QC import app_qc
from genevue.Remote import app_remote
from genevue.Sequences import app_sequence
from genevue.Taxonomy import app_taxonomy
from genevue.Utils import app_file
from genevue.Utils.Completion import install_completion, uninstall_completion

COPYRIGHT = """
Copyright (c) 2026 HYLi360.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

In addition to the above, you must comply with the following
additional terms:

1. Pursuant to Section 7(b), when redistributing this software, you
   must retain the original authors' copyright notice
   (COPYRIGHT in /src/genevue/cli.py), and ensure that all users
   can easily view it.
2. Pursuant to Section 7(c), adapted versions of this software must
   use a version number different from that of the original version
   to distinguish them.
3. Pursuant to Sections 7(a) and 7(f), the software authors make no
   warranty regarding the reliability of this program. Redistributors
   are solely responsible for any warranties they provide.
"""


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
    def apikey(
        provider: Annotated[
            Literal["NCBI", "ENA"],
            typer.Argument(..., help=_("Service provider (e.g. NCBI)")),
        ],
        key: Annotated[str, typer.Argument(..., help=_("API key"))],
    ):
        configure.set_apikey(str(provider), key)
        configure.save()

    @app_config.command(help=_("Set contact email."))
    def email(
        your_email: Annotated[str, typer.Argument(..., help=_("Email address"))],
    ):
        configure.set_email(your_email)
        configure.save()

    @app_config.command()
    def displog_level(
        level: Annotated[
            Literal["DEBUG", "INFO", "WARNING", "ERROR"], typer.Argument(..., help="")
        ],
    ):
        configure.set_displog_level(str(level))
        configure.save()

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

    app.add_typer(app_qc, name="qc")

    pipeline = typer.Typer(help=_(" "))
    app.add_typer(pipeline, name="pipeline")

    app.add_typer(app_sequence, name="sequence")

    app.add_typer(app_taxonomy, name="taxonomy")

    app.add_typer(app_remote, name="remote")

    app.add_typer(app_file, name="file")

    app_config = _build_config_commands(_, configure)
    app.add_typer(app_config, name="config")

    app.add_typer(app_preprocessing, name="preprocessing")

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

    @app.command(name="version", help=_("Print genevue version."))
    def cmd_version():
        d = Diagnosis()

        console.print(f"\n[green]GeneVue[/green]    [cyan]{__full_version__}[/cyan]\n")
        console.print("This program based on:")
        d.export_simple()
        console.print(
            "Copyright (c) 2026 HYLi360. All rights reserved.", highlight=False
        )
        console.print("For details, please execute 'genevue copyright'.\n")

    @app.command(name="copyright", help=_("Print copyright infomation."))
    def cmd_copyright():
        console.print(COPYRIGHT, highlight=False)

    return app
