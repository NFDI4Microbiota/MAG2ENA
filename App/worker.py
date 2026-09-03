"""Start an RQ worker with MAG2ENA's import path configured reliably."""

from __future__ import annotations

import os
import sys
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent
os.chdir(APP_DIR)
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

from redis import Redis  # noqa: E402
from rq import Queue, Worker  # noqa: E402

# Fail during startup instead of after a user job has been accepted.
from modules.modify import modify_task  # noqa: E402, F401
from modules.submit import submission_task  # noqa: E402, F401
from utils.tasks import _on_failure, _on_success  # noqa: E402, F401


def main() -> None:
    connection = Redis(host="localhost", port=6379)
    queue = Queue("ena", connection=connection)
    Worker([queue], connection=connection).work()


if __name__ == "__main__":
    main()
