# Review

## Scope

- Re-reviewed the summary contract gap analysis using local reference package lessons from `repo/CandiHap`, `repo/geneHapR`, and `repo/hastat`.
- No product source code was changed.

## Recommendation Summary

Critical:

- None.

High:

- Promote display haplotype IDs, display-ready allele states, per-haplotype samples, and explicit count/frequency labels into the authoritative summary contract.
- Make population-aware frequency summaries a first-class summary output when `--population` is provided.

Medium:

- Add explicit selector/region metadata to summary output.
- Move per-site functional category into summary only if it is intended as user result information, not merely a plot strip.
- Rename or split `freq`, because the current value is a count, not a frequency.

Low:

- Keep indel aliases plot-side unless the aliases are also intended to appear in user-facing result tables.
- Do not optimize plot style/layout in this pass.

## Reference Rationale

- `geneHapR` treats haplotype summaries, accessions, phenotype/population relationships, and visualization as connected but user-visible domain outputs. This supports making `summary` richer rather than leaving key result meaning hidden in plotting code.
- `CandiHap`'s strongest lesson is product clarity: users need a complete candidate-gene report, not only a figure. This supports a document-like summary table.
- `hastat` shows useful module separation, but also demonstrates the risk of plot/config code becoming a secondary data contract. This supports keeping plot selective and summary authoritative.

## External Model Review

CCG dual-model analysis/review was not run because `~/.claude/bin/codeagent-wrapper` is not present on this machine.
