#!/usr/bin/env python3
"""Create a new session log file."""

from datetime import datetime
from pathlib import Path

def create_session():
    # Use absolute path relative to project root or current working dir
    log_dir = Path(".context/memories/session_logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    
    today = datetime.now().strftime("%Y-%m-%d")
    
    # Find existing sessions for today
    existing = list(log_dir.glob(f"{today}-session-*.md"))
    session_num = len(existing) + 1
    
    filename = f"{today}-session-{session_num:02d}.md"
    filepath = log_dir / filename
    
    template = f"""# Session Log: {today} (Session {session_num})

**Date**: {today}
**Time**: {datetime.now().strftime("%H:%M")} - ...
**Focus**: ...

---

## Key Topics
- ...

---

## Decisions Made
- ...

---

## Action Items
| Action | Owner | Status |
|--------|-------|--------|
| ... | ... | Pending |

---

## Session Closed
**Status**: Open
"""
    
    filepath.write_text(template)
    print(f"✅ Created: {filepath}")
    print(f"   Session: {today}-session-{session_num:02d}")

if __name__ == "__main__":
    create_session()
