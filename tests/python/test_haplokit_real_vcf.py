from __future__ import annotations

import json
import sys
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit.cli import main


DATA_DIR = ROOT / "data"


@pytest.fixture()
def indexed_vcf(tmp_path: Path) -> Path:
    src = DATA_DIR / "var.vcf.gz"
    decompressed = tmp_path / "var.vcf"
    sorted_vcf = tmp_path / "var.sorted.vcf"
    with pysam.BGZFile(str(src), "rb") as r, decompressed.open("wb") as w:
        w.write(r.read())
    header: list[str] = []
    records: list[tuple[str, int, str]] = []
    for line in decompressed.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            header.append(line)
            continue
        chrom, pos, *_ = line.split("\t")
        records.append((chrom, int(pos), line))
    records.sort(key=lambda item: (item[0], item[1]))
    sorted_vcf.write_text(
        "\n".join(header + [record[2] for record in records]) + "\n",
        encoding="utf-8",
    )
    compressed = tmp_path / "var.repacked.vcf.gz"
    pysam.tabix_compress(str(sorted_vcf), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def _read_tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _summary_meta(path: Path) -> dict[str, object]:
    rows = _read_tsv(path)
    site_positions = rows[1][1:-2]
    return {
        "haplotypes": int(rows[0][-1]),
        "individuals": int(rows[1][-1]),
        "variants": int(rows[2][-1]),
        "site_positions": site_positions,
        "summary_cols_tail": rows[3][-2:],
    }


def test_region_mode_writes_hapresult_and_hap_summary_tsv(tmp_path: Path, indexed_vcf: Path) -> None:
    out_dir = tmp_path / "out"
    exit_code = main(["view", str(indexed_vcf), "-r", "scaffold_1:4300-5000", "--output-file", str(out_dir)])

    assert exit_code == 0
    hapresult = out_dir / "hapresult.tsv"
    hap_summary = out_dir / "hap_summary.tsv"
    assert hapresult.exists()
    assert hap_summary.exists()

    meta = _summary_meta(hap_summary)
    assert meta["variants"] > 0
    assert meta["individuals"] > 0
    assert meta["haplotypes"] > 0
    assert meta["summary_cols_tail"] == ["Accession", "freq"]

    result_rows = _read_tsv(hapresult)
    assert result_rows[3][-1] == "Accession"
    assert len(result_rows) > 4


def test_site_mode_with_r_chr_pos_behaves_as_one_site(tmp_path: Path, indexed_vcf: Path) -> None:
    out_dir = tmp_path / "out_site"
    exit_code = main(["view", str(indexed_vcf), "-r", "scaffold_1:4300", "--output-file", str(out_dir)])

    assert exit_code == 0
    meta = _summary_meta(out_dir / "hap_summary.tsv")
    assert meta["variants"] == 1
    assert len(meta["site_positions"]) == 1
    assert meta["site_positions"][0] == "4300"


def test_samples_file_limits_individual_count_in_tsv(tmp_path: Path, indexed_vcf: Path) -> None:
    sample_file = tmp_path / "samples.list"
    sample_file.write_text("C1\nC5\nC16\n", encoding="utf-8")
    out_dir = tmp_path / "out_subset"

    exit_code = main(
        [
            "view",
            str(indexed_vcf),
            "-r",
            "scaffold_1:4300-5000",
            "-S",
            str(sample_file),
            "--output-file",
            str(out_dir),
        ]
    )

    assert exit_code == 0
    meta = _summary_meta(out_dir / "hap_summary.tsv")
    assert meta["individuals"] == 3


def test_output_file_prefix_yields_prefixed_hapresult_and_hap_summary(tmp_path: Path, indexed_vcf: Path) -> None:
    out_prefix = tmp_path / "custom.tsv"
    exit_code = main(["view", str(indexed_vcf), "-r", "scaffold_1:4300-5000", "--output-file", str(out_prefix)])
    assert exit_code == 0

    assert (tmp_path / "custom.hapresult.tsv").exists()
    assert (tmp_path / "custom.hap_summary.tsv").exists()


def test_bed_samples_plot_gff_pipeline_writes_expected_artifacts(tmp_path: Path, indexed_vcf: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("scaffold_1\t4299\t4300\twin1\nscaffold_1\t4344\t4345\twin2\n", encoding="utf-8")
    sample_file = tmp_path / "samples.list"
    sample_file.write_text("C1\nC5\nC16\n", encoding="utf-8")
    out_dir = tmp_path / "out_bed"

    exit_code = main(
        [
            "view",
            str(indexed_vcf),
            "-R",
            str(bed),
            "-S",
            str(sample_file),
            "--plot",
            "--plot-format",
            "pdf",
            "--gff",
            str(DATA_DIR / "annotation.gff"),
            "--output-file",
            str(out_dir),
        ]
    )
    assert exit_code == 0

    summary_files = sorted(out_dir.glob("hap_summary_*.tsv"))
    result_files = sorted(out_dir.glob("hapresult_*.tsv"))
    plot_files = sorted(out_dir.glob("*.pdf"))
    assert len(summary_files) == 2
    assert len(result_files) == 2
    assert len(plot_files) == 2

    for summary_file in summary_files:
        meta = _summary_meta(summary_file)
        assert meta["individuals"] == 3
        assert meta["summary_cols_tail"] == ["Accession", "freq"]

    gff_summary = out_dir / "gff_ann_summary.tsv"
    assert gff_summary.exists()
    gff_lines = gff_summary.read_text(encoding="utf-8").splitlines()
    assert gff_lines[0] == "chr\tstart\tend\tann"
    assert len(gff_lines) == 3


def test_jsonl_mode_remains_available_for_metadata_regression(tmp_path: Path, indexed_vcf: Path) -> None:
    out_file = tmp_path / "result.jsonl"
    exit_code = main(
        [
            "view",
            str(indexed_vcf),
            "-r",
            "scaffold_1:4300-5000",
            "--max-diff",
            "0.2",
            "--output-format",
            "jsonl",
            "--output-file",
            str(out_file),
        ]
    )

    assert exit_code == 0
    row = json.loads(out_file.read_text(encoding="utf-8").strip())
    assert row["grouping_mode"] == "approx-region"
    assert row["grouping_method"] == "max-diff"
