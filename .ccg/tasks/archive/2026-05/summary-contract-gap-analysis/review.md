# Review

## Scope

- Compared C++ `ViewResult` / JSON output with Python CLI TSV generation and Python plotting transformations.
- Identified Python-derived information that may belong in authoritative summary output.
- Saved the decision aid in `summary-contract-gap-analysis.md`.
- No product source code was changed.

## Verification

- `.ccg/spec` is absent, so no local spec files applied.
- Confirmed C++ summary fields in `src/cpp/view_backend.h` and `serialize_view_result_json()`.
- Confirmed TSV summary generation in `haplokit/cli.py`.
- Confirmed plot-derived transformations in `haplokit/_transform.py`, `haplokit/_table.py`, `haplokit/_distribution.py`, and `haplokit/_network.py`.

## External Model Review

CCG dual-model analysis/review was not run because `~/.claude/bin/codeagent-wrapper` is not present on this machine.

## Findings

- Critical: none.
- Warning: this is a contract gap analysis only; implementation choices are pending owner decision.
- Info: plot layout should remain unchanged until summary contract decisions are made.
