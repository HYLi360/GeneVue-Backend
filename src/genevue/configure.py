import json
import os
from pathlib import Path
from typing import Literal

from rich import print

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
    "PATH": {},
    "API_KEY": {"NCBI": ""},
    "E-MAIL": "",
}


class Configure:
    def __init__(self):
        self.config_path: Path = CONFIG_PATH
        if not CONFIG_PATH.exists():
            print("Configure file not found.")
            print(
                f"Generate configure file to {HOME_PATH / 'genevue' / 'config.json'}."
            )
            json.dump(CONFIG_DEFAULT, open(HOME_PATH / "genevue" / "config.json", "w"))
        self.config: dict = json.load(open(HOME_PATH / "genevue" / "config.json"))

    def reset(self):
        json.dump(CONFIG_DEFAULT, open(HOME_PATH / "genevue" / "config.json", "w"))
        self.config: dict = json.load(open(HOME_PATH / "genevue" / "config.json"))

    def set_apikey(self, service_provider: Literal["NCBI"], new_api_key: str):
        self.config["API_KEY"][service_provider] = new_api_key

    def get_apikey(self, service_provider: Literal["NCBI"]):
        return self.config["API_KEY"][service_provider]

    def set_email(self, new_email: str):
        self.config["E-MAIL"] = new_email
        json.dump(self.config, open(HOME_PATH / "genevue" / "config.json", "w"))

    def get_email(self):
        return self.config["E-MAIL"]

    @property
    def email(self):
        return self.config["E-MAIL"]

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
        return self.config["localization"]["translate_path"]
