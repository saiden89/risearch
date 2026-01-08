#!/usr/bin/env python3
"""Append a checkpoint to the current session log."""

import sys
from datetime import datetime
from pathlib import Path


def quicksave(summary: str):
    log_dir = Path(".context/memories/session_logs")
    today = datetime.now().strftime("%Y-%m-%d")

    # Find today's session logs
    logs = sorted(log_dir.glob(f"{today}-session-*.md"))
    if not logs:
        print(f"⚠️ No session log found for {today}")
        return

    current_log = logs[-1]  # Most recent
    timestamp = datetime.now().strftime("%H:%M")

    checkpoint = f"\n\n### ⚡ Checkpoint [{timestamp}]\n{summary}\n"

    with open(current_log, "a") as f:
        f.write(checkpoint)

    print(f"✅ Quicksave [{timestamp}] → {current_log.name}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        quicksave(" ".join(sys.argv[1:]))
    else:
        print("Usage: python quicksave.py <summary>")
