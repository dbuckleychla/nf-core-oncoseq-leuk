"""Count-only stage logging: never print sample IDs, alleles or clinical paths."""
from datetime import datetime, timezone


def progress(message):
    print('[review ' + datetime.now(timezone.utc).isoformat(timespec='seconds') + '] ' + message,
          flush=True)
