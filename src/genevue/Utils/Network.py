import asyncio
import time
from dataclasses import dataclass, field
from pathlib import Path

import aiohttp
import requests
from rich.progress import Progress

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)

# Retry configuration for AsyncDownloadManager
_RETRY_DELAYS = [1.0, 2.0, 4.0, 4.0, 4.0]
_DEFAULT_CONNECT_TIMEOUT = 30
_DEFAULT_READ_TIMEOUT = 60
_DEFAULT_TOTAL_TIMEOUT = 300


class DownloadFailedError(Exception):
    """Raised when a download fails after exhausting all retries."""


def _pretty_status_display(received_byte: int, spend_time_by_a_chunk: float):
    """Human-readable display of download progress (reused by sync & async managers)."""
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


@dataclass
class DownloadManager:

    def __init__(self, url: str, headers: dict, target_file_path: Path):
        self.url = url
        self.headers = headers
        self.target_file_path = target_file_path
        self.no_total_file_size = True
        self.target_file_path_completed = None

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
        r = ""
        s = ""

        if total_size == 0:
            with console.status("Downloading...") as status:
                for chunk in stream.iter_content(chunk_size=2048):
                    received += f_out.write(chunk)
                    r, s = _pretty_status_display(
                        received, received / (time.time() - start_time)
                    )
                    status.update(f"Received: {r}   Total Speed: {s}")
        else:
            with Progress() as progress:
                task_download = progress.add_task(
                    "[green]Downloading...", total=total_size
                )
                for chunk in stream.iter_content(chunk_size=2048):
                    received += f_out.write(chunk)
                    r, s = _pretty_status_display(
                        received, received / (time.time() - start_time)
                    )
                    progress.update(task_download, advance=received)
        logger.info(f"\n{r} received in {time.time() - start_time}s. Speed: {s}")

        f_out.close()


@dataclass
class AsyncDownloadManager:
    """Async downloader with exponential-backoff retry for transient network errors.

    Retries on: ClientPayloadError (ChunkedEncodingError), ServerConnectionError,
    ServerTimeoutError, asyncio.TimeoutError, ClientOSError.
    Does NOT retry on: 4xx (except 429), InvalidURL.
    """

    url: str
    headers: dict
    target_file_path: Path
    retry_delays: list[float] = field(default_factory=lambda: list(_RETRY_DELAYS))
    connect_timeout: float = _DEFAULT_CONNECT_TIMEOUT
    read_timeout: float = _DEFAULT_READ_TIMEOUT
    total_timeout: float = _DEFAULT_TOTAL_TIMEOUT

    @staticmethod
    async def _stream_to_file(response, path: Path) -> Path:
        received = 0
        start_time = time.time()
        with open(path, "wb") as f:
            async for chunk in response.content.iter_chunked(2048):
                if not chunk:
                    continue
                f.write(chunk)
                received += len(chunk)
                r, s = _pretty_status_display(
                    received, received / (time.time() - start_time)
                )
        logger.info(f"Downloaded {r} in {time.time() - start_time:.1f}s to {path}")
        return path

    async def download(self) -> Path:
        resolved = self.target_file_path.resolve()
        max_attempts = len(self.retry_delays) + 1
        last_exception = None

        for attempt in range(max_attempts):
            try:
                timeout = aiohttp.ClientTimeout(
                    total=self.total_timeout,
                    connect=self.connect_timeout,
                    sock_read=self.read_timeout,
                )
                connector = aiohttp.TCPConnector(force_close=True)
                async with aiohttp.ClientSession(
                    headers=self.headers, timeout=timeout, connector=connector
                ) as session:
                    async with session.get(self.url) as response:
                        response.raise_for_status()
                        return await self._stream_to_file(response, resolved)

            except aiohttp.ClientResponseError as e:
                if e.status == 429 and attempt < len(self.retry_delays):
                    delay = self.retry_delays[attempt] * 2
                    logger.warning(f"Rate limited (429). Retrying in {delay}s...")
                    resolved.unlink(missing_ok=True)
                    await asyncio.sleep(delay)
                    last_exception = e
                    continue
                raise

            except (
                aiohttp.ClientPayloadError,
                aiohttp.ServerConnectionError,
                aiohttp.ServerTimeoutError,
                asyncio.TimeoutError,
                aiohttp.ClientOSError,
            ) as e:
                last_exception = e
                if attempt < len(self.retry_delays):
                    delay = self.retry_delays[attempt]
                    logger.warning(
                        f"Download attempt {attempt + 1} failed: {e}. "
                        f"Retrying in {delay}s..."
                    )
                    resolved.unlink(missing_ok=True)
                    await asyncio.sleep(delay)
                else:
                    logger.error(f"Download failed after {max_attempts} attempts: {e}")

        raise DownloadFailedError(
            f"Download from {self.url} failed after {max_attempts} attempts. "
            f"Last error: {last_exception}"
        )
