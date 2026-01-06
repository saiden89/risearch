---
applyTo: '**'
---
---
# MCP Operational Directives

This document defines the specialized protocols for using the Model Context Protocol (MCP) servers available in this environment. Strict adherence to these protocols is mandatory for precision, safety, and efficiency.

## 1. Code Intelligence & Analysis (`narsil-mcp`)

**Primary Use:** Deep system interrogation, architectural discovery, and security auditing.

- **Neural Semantic Search**: Find similar code or concept via `narsil-mcp_search_code` or `semantic_search`.
- **Structural Mapping**: Map code hierarchies using `narsil-mcp_get_project_structure` and `get_chunks`.
- **Call Graph Discovery**: Use `narsil-mcp_get_call_graph` and `get_callers`/`get_callees` to trace execution paths.
- **Security Audit**: Perform `narsil-mcp_check_owasp_top10` and `trace_taint` for vulnerability scanning.
- **Hotspots**: Identify refactoring targets using `narsil-mcp_get_function_hotspots`.
- **Fuzzy Symbol Recovery**: For legacy code (C/C++) or non-idiomatic declarations where formal symbol indexing (`find_symbols`) fails, ALWAYS use `workspace_symbol_search` or `search_code` to ensure full coverage.

## 2. Rust Ecosystem & Toolchain (`cargo-mcp`)

**Primary Use:** Managing the Rust build cycle and dependencies.

- **Build & Check**: Always use `cargo-mcp_cargo_check` and `cargo_build` instead of raw shell commands.
- **Dependency Management**: Use `cargo_add` and `cargo_remove` to modify `Cargo.toml`.
- **Verification**: Run tests via `cargo-mcp_cargo_test`.
- **Absolute**: never ever try to run cargo commands with run_command, always use the carbo MCP.

## 3. External Knowledge & Documentation (`context7`)

**Primary Use:** Just-in-time acquisition of library schemas and best practices.

- **Prioritize Evidence**: If a library's behavior is ambiguous, use `context7_resolve-library-id` followed by `query-docs`.
- **Zero-Guessing**: Never assume an API's signature. Verify with `context7` before implementation.

## 4. Reflexive Problem Solving (`sequential-thinking`)

**Primary Use:** Complexity management and error recovery.

- **System 2 Thinking**: Use `sequential-thinking` for multi-step diagnostic reasoning when a parity test fails or a design decision is complex.
- **Hypothesis Verification**: Explicitly state and test hypotheses using sequential thoughts before acting.

## 5. Investigation & Document Operations (`bash`, `filesystem`, `ripgrep`, `pdf-reader`)

**Primary Use:** Raw I/O, searching, and document extraction.

- **PDF Extraction**: Use `pdf-reader_read_pdf` for any PDF processing.
  - ❌ NEVER attempt to `cat` or `strings` a PDF file in the shell.
  - ✅ Use `include_full_text` for complete reading or `pages` for targeted extraction.
- **Avoid "The Shell Trap"**:
  - ❌ NO `cat` or `grep` for code search. Use `narsil-mcp_search_code` or `mcp_ripgrep_search`.
  - ❌ NO manual file tree parsing. Use `narsil-mcp_get_project_structure`.
  - ❌ NO terminal-based file reading. Use `mcp_filesystem_read_text_file`.
- **Shell Usage**: Reserve `mcp_bash_run` for custom scripts, system calls (e.g., `wc`, `find`), or commands not supported by specialized MCPs.

## 6. General Protocols

1. **Verification-First**: After every action that modifies state, use an MCP tool to verify the outcome (e.g., `ls` or `cargo_check`).
2. **Error Recovery**: If an MCP returns EOF or connection closed, attempt to re-initialize the server or verify the binary path in `mcp_config.json`.
3. **PRAR Alignment**:
    - **Perceive**: Use `narsil-mcp` and `context7`.
    - **Reason**: Use `sequential-thinking`.
    - **Act**: Use `cargo-mcp` and `filesystem`.
    - **Refine**: Use `narsil-mcp` security and hotspot scans.
