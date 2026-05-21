# haplokit

面向 CLI 的单倍型查看工具，提供 bcftools 风格选择器、C++ 后端加速和 Python 绘图。

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

## 安装

```bash
pip install haplokit
```

> 源码构建需要 Linux/WSL、Python 3.10+、C++17 工具链、CMake 3.22+ — 见[贡献开发](#贡献开发)。

## 快速开始

```bash
haplokit view data/var.sorted.vcf.gz -r scaffold_1:4300-5000 --output-file out
```

输出：

- `out/hapresult.tsv` — 逐样本单倍型详情
- `out/hap_summary.tsv` — 单倍型计数汇总

## 使用场景

### 1. 区域查询 — 严格精确分组

识别基因组区域内的所有不同单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --output-file out
```

在 `out/` 中生成 `hapresult.tsv` + `hap_summary.tsv`。每行单倍型展示精确等位基因模式；含杂合或缺失呼叫的样本被排除。

### 2. 单位点查询

分析单个变异位点的单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1450 --output-file out_site
```

`--by` 对 `chr:pos` 选择器自动推断为 `site`。

### 3. 基因注释 + 图表

在单倍型表上叠加基因结构。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --gff genes.gff3 --plot --output-file out
```

`genes.gff3` 格式（标准 GFF3）：

```text
chr1	.	gene	1000	3000	.	+	.	ID=gene1;Name=GeneA
chr1	.	CDS	1200	1500	.	+	0	ID=cds1;Parent=gene1
```

添加 SnpEff 风格功能分类色带（CDS、UTR、exon、intron、intergenic）到变异位点上。输出图片（`out/*.png`）+ `gff_ann_summary.tsv`。

<img src="plottable.png" alt="单倍型汇总表" width="800">

图表组件：

- **标题**：区域 + 重叠基因名（提供 `--gff` 时）
- **功能色带**（仅 `--gff`）：彩色条按功能类别标注每个变异
- **POS / ALLELE 行**：变异位置和替代等位基因
- **单倍型行**（H001、H002、...）：每位点等位基因；空 = 参考
- **群体列**（`--population`）：各群体各单倍型样本数
- **n/N**：单倍型频率
- **图例**（仅 `--gff`）：功能类别颜色
- **Indel 脚注**：多等位基因 Indel 用上标标记标注

### 4. 群体分组

比较不同群体间的单倍型分布。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --plot --output-file out
```

`popgroup.txt`（Tab 分隔：`sample<TAB>population`）：

```text
C1	wild
C2	wild
C13	landrace
```

在表和图中添加群体列。

### 5. 地理分布图

在采样地点上映射单倍型组成。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --geo data/sample_china_geo.txt --plot --output-file out
```

`sample_china_geo.txt` 和 `sample_world_geo.txt` 是 Tab 分隔的坐标示例（`ID<TAB>longitude<TAB>latitude<TAB>Hap`）。`Hap` 列用于独立绘图示例；CLI 地图绘图会从 VCF 结果中推断每个样本的单倍型。

```text
ID	longitude	latitude	Hap
C1	116.40	39.90	H001
C2	116.40	39.90	H002
C3	116.40	39.90	H001
```

<img src="plotmap.png" alt="单倍型地理分布图" width="600">

`data/` 下包含世界地图示例资源：

- `sample_world_geo.txt` 与 `sample_china_geo.txt` 保持相同的 `ID/Hap` 组成，只把坐标替换为全球采样地点。
- `world_countries.shp`、`world_countries.shx` 和 `world_countries.dbf` 提供示例世界地图 shapefile。
- `sample_world_geo_map.png` 是生成好的世界地图示例图。

<img src="data/sample_world_geo_map.png" alt="世界单倍型地理分布图" width="600">

图表组件：

- **饼图**：每位置单倍型组成；大小 ∝ √(样本数)
- **颜色图例**：单倍型颜色键
- **气泡大小图例**：ggplot2 风格分级圆圈，展示样本数刻度
- **底图**：GeoJSON 省界多边形（中国）或随附的世界地图 shapefile 示例

### 6. 单倍型网络 — popart 风格

构建单倍型网络，以 [popart](https://popart.maths.otago.ac.nz/) (Leigh & Bryant, 2015) 的视觉规范呈现。支持三种推断方法：TCS (Clement et al. 2002)、MSN 和 MJN (Bandelt, Forster & Röhl 1999)。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -p popgroup.txt --network --plot --output-file out
haplokit view in.vcf.gz -r chr1:1000-2000 --network --network-method mjn --plot --output-file out
```

图表组件：

- **节点**：每个单倍型一个圆；面积 ∝ √(样本数)
- **饼图扇区**（配合 `-p`）：单倍型的群体组成
- **边**：理想长度正比于突变距离（力导向布局）
- **边上的横线刻度**：每条横线代表一个突变（popart 规范）
- **小黑点**：TCS 推断的中间（祖先）节点

![网络算法对比 — MSN / TCS / MJN](plotnetwork_3algo.png)

### 7. BED 批处理

一次运行处理多个区域。

```bash
haplokit view in.vcf.gz -R regions.bed --output-file out_batch
```

`regions.bed`（≥3 列，Tab 分隔）：

```text
chr1	1000	2000
chr2	5000	6000
```

每行 BED 独立处理。输出文件按区域 slug 加后缀（`_chr1_1000_2000`）。

### 8. 近似分组

在容差范围内聚类相似单倍型。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 --max-diff 0.2 --output-file out
```

`--max-diff`（0–1）：差异 ≤ 20% 位点的单倍型归为一组。分组模式从 `strict-region` 变为 `approx-region`。

### 9. 样本子集 + 填补

限制分析到特定样本；缺失呼叫按参考等位基因填补。

```bash
haplokit view in.vcf.gz -r chr1:1000-2000 -S samples.list --impute --output-file out
```

`samples.list`（每行一个样本 ID）：

```text
C1
C5
C16
```

`--impute` 将缺失 GT 视为 `0/0`，提高样本保留率。

## 输出文件

### `hapresult.tsv` — 逐样本单倍型详情

```text
CHR     scaffold_1  scaffold_1  ...  Haplotypes:  8
POS     4300        4345        ...  Individuals: 37
INFO    .           .           ...  Variants:    5
ALLELE  G/C         T/A,GG      ...  Accession
H001    G           T           ...  C8;C9;C11;C14;C18;C25;C26;C28;C31;C35
```

- **表头行**（CHR/POS/INFO/ALLELE）：跨列变异元数据
- **单倍型行**（H001–HNNN）：每位点等位基因；空 = 参考；携带该单倍型的样本列表

### `hap_summary.tsv` — 单倍型计数汇总

与 `hapresult.tsv` 表头相同，多一列 `freq`（计数/总数）：

```text
H001  G   T   T   GCCTA  T   10
H002  G   T   T   A      T   8
H003  C   T   T   A      T   8
```

### `gff_ann_summary.tsv` — 基因注释（仅 `--gff`）

```text
chr           start  end   ann
scaffold_1    4300   5000  test1G0387
```

### 图片文件（`--plot`）

格式由 `--plot-format` 设定（默认 `png`）。按区域 slug 命名：`<prefix>.<chr>_<start>_<end>.png`。

## 完整参数

```
haplokit view <input_vcf> (-r <region> | -R <regions.bed>) [options]
```

`<input_vcf>` 须为已索引 VCF/BCF（`.vcf.gz` + `.tbi`，或 BCF 索引）。

| 参数 | 类型 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `-r, --region` | 字符串 | — | `chr:start-end` 或 `chr:pos` |
| `-R, --regions-file` | 路径 | — | BED 文件（≥3 列，Tab 分隔） |
| `-S, --samples-file` | 路径 | — | 每行一个样本 ID |
| `--by` | `auto\|region\|site` | `auto` | 分组模式；auto 按选择器形态推断 |
| `--impute` | 开关 | 关 | 缺失 GT 按参考等位基因填补 |
| `-g, --gff` | 路径 | — | GFF3/GTF 基因注释 |
| `-p, --population` | 路径 | — | Tab 分隔的 样本→群体 映射 |
| `--output` | `summary\|detail` | `summary` | 仅 JSONL 模式生效；TSV 总是双表 |
| `--output-format` | `tsv\|jsonl` | `tsv` | 输出格式 |
| `--output-file` | 路径 | — | 输出目录、前缀或 JSONL 文件 |
| `--plot` | 开关 | 关 | 生成单倍型表图 |
| `--plot-format` | `png\|pdf\|svg\|tiff` | `png` | 图片格式 |
| `--max-diff` | 浮点 [0,1] | — | 近似分组阈值 |
| `--geo` | 路径 | — | 样本地理坐标（用于地图） |
| `--network` | 开关 | 关 | 渲染单倍型网络（popart 风格） |
| `--network-method` | `tcs`/`msn`/`mjn` | `tcs` | 网络推断算法 |

选择器规则：`-r` 与 `-R` 互斥且必须提供其一。`--by site` 仅可用于 `-r chr:pos`。

## 后端

C++ 后端（`haplokit_cpp`）处理 VCF 读取和单倍型分组。发现顺序：

1. `HAPLOKIT_CPP_BIN` 环境变量
2. 打包内二进制：`haplokit/_bin/haplokit_cpp`
3. 仓库构建产物：`build-wsl/haplokit_cpp` → `build/haplokit_cpp`
4. 回退：自动运行 `cmake` 构建

供应商依赖库：

- **[htslib](https://github.com/samtools/htslib)** — VCF/BCF 读取，索引随机访问
- **[gffsub](https://github.com/WWz33/gffsub)** — GFF3/GTF 解析，overlap/nearest-gene 查询

### 网络算法

C++ 实现的单倍型网络算法（MSN、TCS、MJN），支持 SIMD 加速：

- **库**：`libhaplokit_network.a`（C++17）
- **算法**：MSN（最小生成网络）、TCS（统计简约）、MJN（中间连接）
- **优化**：AVX2 SIMD Hamming 距离、OpenMP 并行
- **Python 接口**：`haplokit.network`，自动 C++/Python 回退
- **可视化**：PopART 风格渲染，饼图节点、突变刻度线、性状图例

纯 Python 参考实现归档于 `archive/python_reference_implementation/`。

## 贡献开发

```bash
cmake -S . -B build-wsl && cmake --build build-wsl -j12
HAPLOKIT_CPP_BIN=$PWD/build-wsl/haplokit_cpp python -m pytest -q tests/python
ctest --test-dir build-wsl --output-on-failure
```

## 致谢

设计灵感来自 geneHapR：

> Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9

网络可视化遵循 [popart](https://popart.maths.otago.ac.nz/) 的规范：

> Leigh, J. W. & Bryant, D. popart: full‐feature software for haplotype network construction. Methods in Ecology and Evolution 6, 1110–1116 (2015). https://doi.org/10.1111/2041-210X.12410

## 许可证

GPL-3.0-or-later
