import json
import os
from pathlib import Path
from typing import Literal, Optional

from Bio.SeqUtils.ProtParamData import pa

HOME_PATH = Path(
    os.environ["APPDATA"] if os.name == "nt" else os.path.expanduser("~/.config")
)
Path(HOME_PATH / "genevue").mkdir(exist_ok=True)
CONFIG_PATH = Path(HOME_PATH / "genevue" / "config.json")

TRANSLATE_PATH = (
    Path(os.environ["APPDATA"]) / "genevue" / "locales"
    if os.name == "nt"
    else os.path.expanduser("~/.local/share/genevue/locales")
)

CONFIG_DEFAULT = {
    "localization": {
        "language": "zh_CN",
        "translate_path": TRANSLATE_PATH,
    },
    "LOGLEVEL": 20,
    "PATH": {},
    "API_KEY": {"NCBI": ""},
    "E-MAIL": "",
}


class Configure:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True
        self.config_path: Path = CONFIG_PATH
        if not CONFIG_PATH.exists():
            print("Configure file not found.")
            print(
                f"Generate configure file to {HOME_PATH / 'genevue' / 'config.json'}."
            )
            self._write(CONFIG_DEFAULT)
        self.config = json.load(open(self.config_path))

    def _write(self, config_content: dict) -> None:
        self.config_path.write_text(json.dumps(config_content, indent=2))

    def save(self) -> None:
        """Persist current in-memory config to disk."""
        self._write(self.config)

    def reset(self) -> None:
        self._write(CONFIG_DEFAULT)
        self.config: dict = CONFIG_DEFAULT

    def set_apikey(self, provider: str, new_api_key: str) -> None:
        self.config["API_KEY"][provider] = new_api_key

    def get_apikey(self, provider: Literal["NCBI"]) -> Optional[str]:
        return self.config["API_KEY"].get(provider, None)

    def set_program_path(self, program_name: str, program_path: str):
        self.config["PATH"][program_name] = program_path

    def get_program_path(self, program_name: str) -> Path:
        return Path(self.config["PATH"].get(program_name, program_name))

    @property
    def email(self) -> Optional[str]:
        return self.config.get("E-MAIL", None)

    def set_email(self, new_email: str) -> None:
        self.config["E-MAIL"] = new_email

    @property
    def displog_level(self) -> str:
        return self.config.get("LOGLEVEL", "INFO")

    def set_displog_level(self, new_level: str) -> None:
        self.config["LOGLEVEL"] = new_level

    def set_l10n(self, mode: Literal["change_lang", "change_ts_path"], new_value: str):
        match mode:
            case "change_lang":
                self.config["localization"]["language"] = new_value
            case "change_ts_path":
                self.config["localization"]["translate_path"] = new_value

    @property
    def language(self):
        return self.config["localization"]["language"]

    @property
    def ts_path(self):
        return Path(self.config["localization"]["translate_path"])
