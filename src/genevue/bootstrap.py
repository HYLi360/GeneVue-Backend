from __future__ import annotations

from rich.console import Console
from rich.traceback import install as install_rich_traceback

from genevue.configure import Configure
from genevue.l10n import init_translations
from genevue.logsystem import _setup_busco_bridge


def bootstrap() -> tuple:
    configure = Configure()

    _ = init_translations(
        lang=configure.language,
        ts_path=configure.ts_path,
    )

    console = Console()
    install_rich_traceback(show_locals=True)

    _setup_busco_bridge(console)

    from genevue.cli import _build_app

    app = _build_app(_, configure)

    return app, _, configure, console
