# haplokit

CLI haplotype viewer with bcftools-like selectors, C++ backend, and Python plotting.

<!-- README-I18N:START -->

**English** | [汉语](./README.zh-CN.md)

<!-- README-I18N:END -->

## Installation

```bash
pip install haplokit
```

> Source build requires Linux/WSL, Python 3.10+, C++17 toolchain, CMake 3.22+ — see [Contributing](#contributing).

## Quick Start

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

Output:

- `out/hapresult.tsv` — per-sample haplotype detail
- `out/hap_summary.tsv` — haplotype count summary

## Usage Scenarios

### 1. Region query — strict haplotype grouping

Identify all distinct haplotypes in a genomic region.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
```

Produces `hapresult.tsv` + `hap_summary.tsv` in `out/`. Each haplotype row shows the exact allele pattern; samples with any heterozygous or missing call are excluded.

### 2. Single-site query

Analyze haplotype at one variant position.

```bash
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

`--by` auto-resolves to `site` for `chr:pos` selectors.

### 3. Gene annotation + figure

Overlay gene structure on the haplotype table.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
```

`genes.gff3` format (standard GFF3):

```text
chr1	.	gene	1000	3000	.	+	.	ID=gene1;Name=GeneA
chr1	.	CDS	1200	1500	.	+	0	ID=cds1;Parent=gene1
```

Adds SnpEff-style functional category strip (CDS, UTR, exon, intron, intergenic) above variant positions. Writes figure (`out/*.png`) + `gff_ann_summary.tsv`.

<img src="plottable.png" alt="Haplotype summary table" width="800">

Figure components:

- **Title**: region + overlapping gene name (when `--gff` provided)
- **Function strip** (`--gff` only): colored bar classifying each variant by functional category
- **POS / ALLELE rows**: variant positions and alternate alleles
- **Haplotype rows** (H001, H002, ...): allele per position; empty = reference
- **Population columns** (`--population`): sample counts per haplotype per group
- **n/N**: haplotype frequency
- **Legend** (`--gff` only): functional category colors
- **Indel footnotes**: multi-allele indels annotated with superscript markers

### 4. Population grouping

Compare haplotype distributions across populations.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

`popgroup.txt` (tab-separated: `sample<TAB>population`):

```text
C1	wild
C2	wild
C13	landrace
```

Adds population columns to the table and figure.

### 5. Geographic distribution map

Map haplotype composition at sampling locations.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
```

`sample_china_geo.txt` and `sample_world_geo.txt` are tab-separated coordinate examples (`ID<TAB>longitude<TAB>latitude<TAB>Hap`). The `Hap` column is included for standalone plotting examples; CLI map plotting derives each sample's haplotype from the VCF result.

```text
ID	longitude	latitude	Hap
C1	116.40	39.90	H001
C2	116.40	39.90	H002
C3	116.40	39.90	H001
```

<img src="plotmap.png" alt="Haplotype geographic distribution" width="600">

World example resources are included under `data/`:

- `sample_world_geo.txt` keeps the same `ID/Hap` composition as `sample_china_geo.txt`, but replaces coordinates with global sampling locations.
- `world_countries.shp`, `world_countries.shx`, and `world_countries.dbf` provide the example world map shapefile.
- `sample_world_geo_map.png` is the generated world map example.

<img src="data/sample_world_geo_map.png" alt="World haplotype geographic distribution" width="600">

Figure components:

- **Pie charts**: haplotype composition per location; size ∝ √(sample count)
- **Color legend**: haplotype color key
- **Bubble-size legend**: ggplot2-style graduated circles, showing the sample-count scale
- **Base map**: GeoJSON province boundaries (China) or the bundled world shapefile example

### 6. Haplotype network — popart-style

Build a haplotype network and visualize it in the conventions of [popart](https://popart.maths.otago.ac.nz/) (Leigh & Bryant 2015). Supports three inference methods: TCS (Clement et al. 2002), MSN and MJN (Bandelt, Forster & Röhl 1999).

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out
haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
```

Figure components:

- **Nodes**: one circle per haplotype; area ∝ √(sample count)
- **Pie slices** (with `-p`): population composition per haplotype
- **Edges**: ideal length proportional to mutation distance (force-directed layout)
- **Hatch marks across edges**: one tick per mutation (popart convention)
- **Small black dots**: inferred median (intermediate) vertices, where TCS infers ancestors

![Network algorithms comparison — MSN / TCS / MJN](plotnetwork_3algo.png)

### 7. BED batch processing

Process multiple regions in one run.

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

`regions.bed` (≥3 tab-separated columns):

```text
chr1	1000	2000
chr2	5000	6000
```

Each BED row is processed independently. Output files are suffixed by region slug (`_chr1_1000_2000`).

### 8. Approximate grouping

Cluster similar haplotypes within a tolerance.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff` (0–1): haplotypes differing at ≤ 20% of positions merge into one group. Grouping mode changes from `strict-region` to `approx-region`.

### 9. Sample subset + imputation

Restrict analysis to specific samples; fill missing calls as reference.

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list` (one sample ID per line):

```text
C1
C5
C16
```

`--impute` treats missing GT as `0/0`, increasing sample retention.

## Output Files

### `hapresult.tsv` — per-sample haplotype detail

```text
CHR     scaffold_1  scaffold_1  ...  Haplotypes:  8
POS     4300        4345        ...  Individuals: 37
INFO    .           .           ...  Variants:    5
ALLELE  G/C         T/A,GG      ...  Accession
H001    G           T           ...  C8;C9;C11;C14;C18;C25;C26;C28;C31;C35
```

- **Header rows** (CHR/POS/INFO/ALLELE): variant metadata across columns
- **Haplotype rows** (H001–HNNN): allele at each position; empty = reference; list of samples carrying this haplotype

### `hap_summary.tsv` — haplotype count summary

Same header as `hapresult.tsv`, plus a `freq` column (count/total):

```text
H001  G   T   T   GCCTA  T   10
H002  G   T   T   A      T   8
H003  C   T   T   A      T   8
```

### `gff_ann_summary.tsv` — gene annotation (`--gff` only)

```text
chr           start  end   ann
scaffold_1    4300   5000  test1G0387
```

### Figure files (`--plot`)

Format set by `--plot-format` (default `png`). Named per region slug: `<prefix>.<chr>_<start>_<end>.png`.

## Full Parameters

```
haplokit view <input_vcf> (-r <region> | -R <regions.bed>) [options]
```

`<input_vcf>` must be an indexed VCF/BCF (`.vcf.gz` + `.tbi`, or BCF index).

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-r, --region` | string | — | `chr:start-end` or `chr:pos` |
| `-R, --regions-file` | path | — | BED file (≥3 tab-separated columns) |
| `-S, --samples-file` | path | — | One sample ID per line |
| `--by` | `auto\|region\|site` | `auto` | Grouping mode; auto infers from selector shape |
| `--impute` | flag | off | Impute missing GT as reference |
| `-g, --gff` | path | — | GFF3/GTF for gene annotation |
| `-p, --population` | path | — | Tab-separated sample → population map |
| `--output` | `summary\|detail` | `summary` | JSONL mode only; TSV always writes both |
| `--output-format` | `tsv\|jsonl` | `tsv` | Output format |
| `--output-file` | path | — | Output directory, prefix, or JSONL file |
| `--plot` | flag | off | Generate haplotype table figure |
| `--plot-format` | `png\|pdf\|svg\|tiff` | `png` | Figure format |
| `--max-diff` | float [0,1] | — | Approximate grouping threshold |
| `--geo` | path | — | Sample geographic coordinates for map |
| `--network` | flag | off | Render haplotype network (popart-style) |
| `--network-method` | `tcs`/`msn`/`mjn` | `tcs` | Network inference algorithm |

Selector rules: `-r` and `-R` are mutually exclusive and one is required. `--by site` only valid with `-r chr:pos`.

## Backend

C++ backend (`haplokit_cpp`) handles VCF reading and haplotype grouping. Discovery order:

1. `HAPLOKIT_CPP_BIN` env var
2. Packaged binary: `haplokit/_bin/haplokit_cpp`
3. Repo build: `build-wsl/haplokit_cpp` → `build/haplokit_cpp`
4. Fallback: auto-run `cmake` build

Vendored libraries:

- **[htslib](https://github.com/samtools/htslib)** — VCF/BCF reading with indexed random access
- **[gffsub](https://github.com/WWz33/gffsub)** — GFF3/GTF parsing with overlap/nearest-gene queries

### Network Algorithms

C++ implementation of haplotype network algorithms (MSN, TCS, MJN) with SIMD acceleration:

- **Library**: `libhaplokit_network.a` (1.7 MB, C++17)
- **Algorithms**: MSN (Minimum Spanning Network), TCS (Statistical Parsimony), MJN (Median-Joining)
- **Optimizations**: AVX2 SIMD Hamming distance, OpenMP parallelization, O(1) edge deletion
- **Status**: Core C++ library compiled and tested ✓
- **Python Interface**: `haplokit.network` with automatic C++/Python fallback
- **Visualization**: PopART-style rendering with pie chart nodes, hatch marks, trait legends

Reference implementation (pure Python) archived in `archive/python_reference_implementation/` for algorithm verification.

Performance (100 haplotypes, 1000bp):
- MSN: ~10ms (C++) vs ~1s (Python)
- TCS: ~20ms (C++) vs ~2s (Python)
- MJN: ~50ms (C++) vs ~5s (Python)

## Contributing

```bash
cmake -S . -B build-wsl && cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

## Acknowledgements

Inspired by geneHapR:

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

Network visualization follows the conventions of [popart](https://popart.maths.otago.ac.nz/):

> Leigh, J. W. & Bryant, D. popart: full‐feature software for haplotype network construction. Methods in Ecology and Evolution 6, 1110–1116 (2015). https://doi.org/10.1111/2041-210X.12410

## License

GPL-3.0-or-later
