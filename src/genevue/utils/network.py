import time

import requests
from pathlib import Path
from dataclasses import dataclass

from genevue import console, setup_rich_logger
from rich.progress import Progress

logger = setup_rich_logger(__name__)


@dataclass
class DownloadManager:

    def __init__(self, url: str, headers: dict, target_file_path: Path):
        self.url = url
        self.headers = headers
        self.target_file_path = target_file_path
        self.no_total_file_size = True
        self.target_file_path_completed = None

    @staticmethod
    def _pretty_status_display(received_byte: int, spend_time_by_a_chunk: float):
        # B, KiB, MiB, GiB, TiB
        received_size_human_ls = [
            (received_byte, "B"),
            (received_byte / 1024, "KiB"),
            (received_byte / 1024 / 1024, "MiB"),
            (received_byte / 1024 / 1024 / 1024, "GiB"),
            (received_byte / 1024 / 1024 / 1024 / 1024, "TiB"),
        ]
        received_size_human = "???"
        for r in received_size_human_ls:
            if r[0] < 1:
                break
            received_size_human = f"{r[0]:.2f} {r[1]}"

        spend_time_by_a_chunk_human_ls = [
            (spend_time_by_a_chunk, "B/s"),
            (spend_time_by_a_chunk / 1024, "KiB/s"),
            (spend_time_by_a_chunk / 1024 / 1024, "MiB/s"),
            (spend_time_by_a_chunk / 1024 / 1024 / 1024, "GiB/s"),
            (spend_time_by_a_chunk / 1024 / 1024 / 1024 / 1024, "TiB/s"),
        ]
        spend_time_by_a_chunk_human = "???"
        for s in spend_time_by_a_chunk_human_ls:
            if s[0] < 1:
                break
            spend_time_by_a_chunk_human = f"{s[0]:.2f} {s[1]}"

        return received_size_human, spend_time_by_a_chunk_human

    def download(self):
        self.target_file_path_completed = self.target_file_path.resolve()
        logger.info(
            f"Start downloading from {self.url} to {self.target_file_path_completed}."
        )
        stream = requests.get(url=self.url, headers=self.headers, stream=True)
        if stream.headers.get("Content-Length") is None:
            logger.info(
                "Unable to determine the file size. Progress bar has been disabled."
            )
            total_size = 0
        else:
            total_size = float(stream.headers.get("Content-Length"))

        f_out = open(self.target_file_path_completed, "wb")
        received = 0
        start_time = time.time()

        if total_size == 0:
            # Status method
            with console.status("Start Downloading.") as status:
                for chunk in stream.iter_content(chunk_size=8192):
                    received += f_out.write(chunk)
                    r, s = self._pretty_status_display(
                        received, received / (time.time() - start_time)
                    )
                    status.update(f"Received: {r}   Total Speed: {s}")

        else:
            with Progress() as progress:
                task_download = progress.add_task(
                    "[green]Downloading...", total=total_size
                )
                for chunk in stream.iter_content(chunk_size=8192):
                    received += f_out.write(chunk)
                    r, s = self._pretty_status_display(
                        received, received / (time.time() - start_time)
                    )
                    progress.update(task_download, advance=received)
        logger.info(f"\n{r} received in {time.time() - start_time}s.")

        f_out.close()
