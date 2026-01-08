---
description: Boot the AI assistant with context
---

# /start — Execution Script

## Phase 1: Load Identity

- [ ] Read `.framework/modules/Core_Identity.md`

## Phase 2: Recall Context

- [ ] Find the latest session log in `.context/memories/session_logs/`
- [ ] Display a summary of the last session

## Phase 3: Create New Session

// turbo

- [ ] Run `python3 .agent/scripts/create_session.py`

## Phase 4: Confirm Ready

- [ ] Output: "⚡ Ready. (Session XX started.)"
