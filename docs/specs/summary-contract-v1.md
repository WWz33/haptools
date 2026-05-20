# Summary Contract v1

`haplokit view --output-format jsonl --output summary` emits the authoritative
user-facing result contract. TSV files and plots are projections of this
summary. They may select fewer fields, but should not invent haplotype labels,
allele display states, sample membership, or frequency meaning independently.

## Scope

This contract covers one selector result row. A selector can be a region, site,
or one row from a BED file.

Stable top-level fields:

| Field | Type | Meaning |
| --- | --- | --- |
| `grouping_mode` | string | Strict or approximate grouping mode, e.g. `strict-region` |
| `grouping_method` | string | Haplotype grouping method, e.g. `exact` or `max-diff` |
| `haplotype_label` | object | Label configuration used for display IDs |
| `selector` | object | Parsed selector payload |
| `region_label` | string | Original display region for this result |
| `sites` | array | Variant display metadata used by table/plot projections |
| `haplotypes` | array | Summary rows, present in summary output |
| `sample_count` | integer | Total samples considered before haplotype exclusion |
| `haplotype_count` | integer | Number of haplotype summary rows |
| `variant_count` | integer | Number of variant sites in the selector |

## Haplotype Labels

Default label configuration:

```json
{"prefix": "Hap", "pad": 2}
```

Default display IDs are `Hap01`, `Hap02`, and so on. Users can override this
with optional CLI flags:

```bash
--hap-prefix H --hap-pad 3
```

That produces IDs such as `H001`.

## Haplotype Rows

Each item in `haplotypes` has these stable fields:

| Field | Type | Meaning |
| --- | --- | --- |
| `id` | string | Display haplotype ID, e.g. `Hap01` |
| `hap` | string | Backward-compatible raw haplotype pattern |
| `pattern` | string | Raw haplotype pattern; same value as `hap` in v1 |
| `states` | array[string] | Display-ready allele states for table/plot output |
| `samples` | array[string] | Samples assigned to this haplotype |
| `count` | integer | Number of samples in this haplotype |
| `total` | integer | Denominator used for `frequency` |
| `frequency` | number | `count / total` |
| `frequency_label` | string | Human-readable ratio, e.g. `10/37` |

When a population file is provided, Python may enrich each haplotype with:

| Field | Type | Meaning |
| --- | --- | --- |
| `populations` | array | Per-population count/frequency summary |

Population items contain `population`, `count`, `total`, `frequency`, and
`frequency_label`.

## Projection Rules

- `hap_summary.tsv` should use `haplotypes[].id`, `states`, `samples`, and
  `count`.
- `hapresult.tsv` should map each accession's raw `hap` through the summary
  label map and use summary `states` where available.
- Plotting should consume summary/detail data and select fields for display.
  Plotting code should not generate core result semantics.
- The legacy `hap` field remains for compatibility. New code should prefer
  `pattern` when it needs the raw haplotype pattern.

## Non-Goals

This contract does not define visual style, layout, or publication figure
composition. Those are separate plotting concerns built on top of the summary.
