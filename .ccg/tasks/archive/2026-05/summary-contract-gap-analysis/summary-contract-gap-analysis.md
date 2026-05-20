# Summary Contract Gap Analysis

Goal: make `summary` the authoritative, user-facing result layer. Plotting can keep its current behavior, but plot-only data enrichment should be identified so the owner can decide what belongs in C++ summary output.

Assumptions:

- `plottable`/current plot table content represents the intended user-facing documentation content.
- Plot rendering style and layout are not part of this task.
- This document is a decision aid only; no C++ or Python product code was changed.

## Current Contract Layers

| Layer | Current source | Stability | Notes |
| --- | --- | --- | --- |
| C++ view JSON | `src/cpp/view_backend.h`, `serialize_view_result_json()` | Medium | Central structs exist, but no formal schema file. |
| CLI TSV summary | `haplokit/cli.py::_write_selector_summary_txt()` | Medium | More user-facing than JSON; generated manually in Python. |
| Plot table input | `haplokit/_transform.py`, `haplokit/_table.py` | Low/implicit | Depends on row labels and columns such as `CHR`, `POS`, `ALLELE`, `Accession`, `freq`. |
| Network JSON | `src/cpp/network/network_backend.cpp`, `haplokit/network.py` | Medium | Simple and tested, but separate from view summary. |

## Current C++ View JSON

Fields emitted by C++:

- `grouping_method`
- `grouping_mode`
- `haplotype_count`
- `imputed_ref`
- `max_diff`
- `output_mode`
- `sample_count`
- `variant_count`
- `sites[]`: `chrom`, `pos`, `allele`
- `haplotypes[]`: `hap`, `count`
- `accessions[]`: `hap`, `sample` only in detail/both mode
- `annotation`: `mode`, and when available `id`, `seqid`, `start`, `end`, `strand`, `distance`

Important limitation: `haplotypes[]` contains raw hap strings and counts only. Display IDs, sample lists, population summaries, and display-ready allele cells are not first-class C++ summary fields.

## Current TSV Summary

`haplokit/cli.py` converts C++ JSON into a more user-facing table:

```text
CHR     <chrom...>      Haplotypes:   <n>
POS     <pos...>        Individuals:  <n>
INFO    <annotation...> Variants:     <n>
ALLELE  <alleles...>    Accession     freq
H001    <states...>     sample1;...   <count>
```

This means the TSV already has data not directly present in C++ summary mode:

- stable-looking display hap IDs: `H001`, `H002`, ...
- display allele states converted from encoded hap strings
- accession/sample list per haplotype
- annotation text in the `INFO` row

These are generated in Python from C++ `both` output.

## Python Plot-Derived Information

### 1. Display haplotype IDs

Current Python source:

- `cli.py::_build_hap_label_map()`
- `cli.py::_write_selector_summary_txt()`
- `cli.py::_write_plot_artifacts()`

Current behavior:

- C++ emits raw hap strings.
- Python assigns `H001`, `H002`, ... based on sorted haplotype order.

Decision point:

- If `H001` labels are part of the user-visible result, they should become explicit summary fields, not only positional Python labels.

Candidate C++ summary field:

```text
haplotypes[].id = "H001"
haplotypes[].pattern = raw encoded hap string
```

### 2. Display-ready allele cells

Current Python source:

- `cli.py::_state_to_label()`
- `cli.py::_hap_states()`

Current behavior:

- C++ emits `hap` as encoded genotype state string.
- Python maps states to display labels using site `allele` values.

Decision point:

- If the user should read allele cells directly from summary, C++ summary should include display allele values or an explicit per-site hap state table.

Candidate C++ summary fields:

```text
haplotypes[].states[] = display-ready allele values
haplotypes[].encoded_states[] = raw genotype state values
```

### 3. Sample list per haplotype

Current Python source:

- C++ emits `accessions[]` in detail/both mode.
- `cli.py::_write_selector_summary_txt()` groups samples by hap and writes them into `Accession`.

Current behavior:

- Raw C++ summary mode lacks per-haplotype sample lists.
- TSV summary includes them only because CLI requests `both` internally.

Decision point:

- If summary is the authoritative result layer, per-haplotype sample lists likely belong in C++ summary output.

Candidate C++ summary field:

```text
haplotypes[].samples[] = sample IDs carrying the haplotype
```

### 4. Frequency count and `n/N`

Current Python source:

- C++ emits `haplotypes[].count`.
- `hap_summary.tsv` writes count as `freq`.
- `_transform.py::transform_for_display()` converts sample list into `n/N`.

Current behavior:

- User-facing plot table shows `n/N`.
- Raw TSV column is named `freq` but contains a count, not a proportion.

Decision point:

- Decide whether summary should expose count, fraction, or both. Current `freq` naming is ambiguous.

Candidate C++ summary fields:

```text
haplotypes[].count = 10
haplotypes[].total = 37
haplotypes[].frequency = 0.270270...
haplotypes[].frequency_label = "10/37"
```

### 5. Population breakdown

Current Python source:

- `_transform.py::read_popgroup()`
- `_transform.py::transform_for_display()`
- `cli.py::_write_plot_artifacts()` for network/map inputs

Current behavior:

- Population file is read in Python.
- Plot table adds one column per population with values like `pc/pop_total`.
- Network plot builds `pop_data_net` from detail accessions and pop group mapping.
- Map plot builds `sample_hap` and geospatial sample records in Python.

Decision point:

- If population-aware summary is expected, summary should include per-haplotype population counts/fractions when `--population` is provided.

Candidate C++/summary fields:

```text
haplotypes[].populations[] = {
  population,
  count,
  total,
  frequency,
  frequency_label
}
```

Open design point:

- C++ currently does not receive the population file. Either pass population mapping into C++, or keep population aggregation in Python but write it into the authoritative summary contract before plotting.

### 6. Region title / selector span

Current Python source:

- `_transform.py::transform_for_display()`
- `cli.py::_selector_span()`
- `cli.py::_compose_row()`

Current behavior:

- JSONL rows include `selector` from Python.
- C++ view result itself does not own a `selector` or `region` object.
- Plot title derives `chrom:min_pos-max_pos` from TSV rows.

Decision point:

- If summary should document exactly what was queried, include selector/region explicitly.

Candidate summary fields:

```text
selector = {type, chrom, start, end, record_id?}
region_label = "scaffold_1:4300-5000"
```

### 7. Annotation display text

Current Python source:

- `cli.py::_annotation_text()`
- `cli.py::_info_cells()`
- `_table.py` re-parses GFF for per-position functional strip

Current behavior:

- C++ annotation object contains overlap/nearest metadata.
- TSV `INFO` row only gets one display text cell.
- Plot functional strip re-parses GFF in Python and classifies each variant position.

Decision point:

- There are two annotation concepts:
  1. region-level nearest/overlap gene summary
  2. per-site functional category used in the plot strip

Candidate summary fields:

```text
annotation.region = {mode, id, seqid, start, end, strand, distance}
sites[].function = CDS | UTR | exon | intron | intergenic
sites[].overlapping_features[] = ...
```

### 8. Indel footnote mapping

Current Python source:

- `_transform.py::_indel_footnotes()`

Current behavior:

- Plot shortens long indel allele strings to `i1`, `i2`, ...
- Footnote is display-only today.

Decision point:

- If shortened indel notation is part of user documentation, summary should define the alias mapping.
- If it is only a plot layout trick, keep it in Python.

Candidate summary field:

```text
indel_aliases[] = {alias, allele}
```

## Python Computation That Should Probably Stay Plot-Side

These are rendering/layout computations, not authoritative biological summary fields:

- color palette selection
- luminance-aware label colors
- table cell dimensions
- map symbol radius scaling
- network node radii for plotting
- spring layout coordinates unless user wants reproducible layout export
- hatch mark positions
- figure export format/DPI decisions

## Python Computation That May Need Promotion

These are currently Python-derived but look like user result information:

| Python-derived information | Current owner | Candidate destination |
| --- | --- | --- |
| `H001` display IDs | Python CLI | C++/summary contract |
| display allele cells | Python CLI | C++/summary contract |
| sample list per haplotype | Python CLI from C++ detail | C++/summary contract |
| `n/N` frequency labels | Python transform | summary contract |
| population count/fraction per haplotype | Python transform/CLI | summary contract |
| queried selector/region label | Python CLI | summary contract |
| per-site functional category | Python plot/GFF parser | summary contract if intended as result |
| indel alias table | Python transform | summary contract only if user-facing |

## Suggested Decision Order

No implementation should begin until these are decided:

1. Should C++ summary `haplotypes[]` include `id`, `states`, and `samples`?
2. Should `freq` remain a count, or should summary expose `count`, `frequency`, and `n/N` separately?
3. Should population aggregation become part of summary when `--population` is provided?
4. Should per-site functional categories move from plot-only GFF parsing into summary?
5. Should indel aliases be a formal output field or remain plot-only?

## Current Recommendation Boundary

Do not change plot layout yet. The first implementation step, if approved, should be contract-first: define the expected summary JSON/TSV fields, then update C++/CLI to emit them, then leave plot behavior as-is or adapt it minimally to consume the richer summary.
