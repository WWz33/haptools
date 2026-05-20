from __future__ import annotations

import json
import sys
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haplokit.cli import build_parser, main


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


def test_view_requires_one_selector_source() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["view"])


def test_view_rejects_r_and_r_file_together() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["view", "-r", "chr1:1-10", "-R", "regions.bed"])


def test_view_defaults_to_tsv_and_infers_region_mode_for_interval() -> None:
    parser = build_parser()
    args = parser.parse_args(["view", "-r", "chr1:1-10"])

    assert args.output_mode == "summary"
    assert args.output_format == "tsv"
    assert args.by == "region"


def test_view_infers_site_mode_from_single_position_region() -> None:
    parser = build_parser()
    args = parser.parse_args(["view", "-r", "chr1:1450"])

    assert args.by == "site"
    assert args.region == "chr1:1450"


def test_view_rejects_conflicting_explicit_by_values() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["view", "-r", "chr1:1450", "--by", "region"])
    with pytest.raises(SystemExit):
        parser.parse_args(["view", "-r", "chr1:1-10", "--by", "site"])


def test_view_accepts_gff_alias() -> None:
    parser = build_parser()
    args = parser.parse_args(["view", "-r", "chr1:1-10", "--gff", "anno.gff"])

    assert args.gff3 == "anno.gff"


def test_main_default_tsv_writes_hapresult_and_hap_summary(tmp_path: Path, indexed_vcf: Path) -> None:
    out_dir = tmp_path / "out"
    exit_code = main(["view", str(indexed_vcf), "-r", "scaffold_1:4300-5000", "--output-file", str(out_dir)])

    assert exit_code == 0
    hapresult = out_dir / "hapresult.tsv"
    hap_summary = out_dir / "hap_summary.tsv"
    assert hapresult.exists()
    assert hap_summary.exists()

    summary_rows = _read_tsv(hap_summary)
    assert summary_rows[0][0] == "CHR"
    assert summary_rows[1][0] == "POS"
    assert summary_rows[2][0] == "INFO"
    assert summary_rows[3][0] == "ALLELE"
    assert summary_rows[3][-2:] == ["Accession", "freq"]
    assert summary_rows[4][0] == "Hap01"

    result_rows = _read_tsv(hapresult)
    assert result_rows[0][0] == "CHR"
    assert result_rows[3][-1] == "Accession"
    assert result_rows[4][0].startswith("Hap")
    assert "Hap01" in {row[0] for row in result_rows[4:]}


def test_main_bed_writes_selector_scoped_tsv_names(tmp_path: Path, indexed_vcf: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("scaffold_1\t4299\t4300\twin1\nscaffold_1\t4344\t4345\twin2\n", encoding="utf-8")
    out_dir = tmp_path / "out"

    exit_code = main(["view", str(indexed_vcf), "-R", str(bed), "--output-file", str(out_dir)])
    assert exit_code == 0

    expected = [
        out_dir / "hapresult_scaffold_1_4299_4300.tsv",
        out_dir / "hapresult_scaffold_1_4344_4345.tsv",
        out_dir / "hap_summary_scaffold_1_4299_4300.tsv",
        out_dir / "hap_summary_scaffold_1_4344_4345.tsv",
    ]
    for item in expected:
        assert item.exists(), f"expected file missing: {item}"


def test_main_still_supports_jsonl_when_explicitly_requested(tmp_path: Path, indexed_vcf: Path) -> None:
    out_file = tmp_path / "result.jsonl"
    exit_code = main(
        [
            "view",
            str(indexed_vcf),
            "-r",
            "scaffold_1:4300-5000",
            "--output-format",
            "jsonl",
            "--output-file",
            str(out_file),
        ]
    )
    assert exit_code == 0
    payload = json.loads(out_file.read_text(encoding="utf-8").strip())
    assert payload["grouping_mode"] == "strict-region"
    assert payload["haplotype_label"] == {"prefix": "Hap", "pad": 2}
    assert payload["region_label"] == "scaffold_1:4300-5000"

    first_hap = payload["haplotypes"][0]
    assert first_hap["id"] == "Hap01"
    assert first_hap["pattern"] == first_hap["hap"]
    assert first_hap["count"] > 0
    assert first_hap["total"] == payload["sample_count"]
    assert first_hap["frequency"] > 0
    assert first_hap["frequency_label"] == f"{first_hap['count']}/{payload['sample_count']}"
    assert first_hap["states"]
    assert first_hap["samples"]


def test_main_jsonl_accepts_optional_haplotype_label_override(tmp_path: Path, indexed_vcf: Path) -> None:
    out_file = tmp_path / "result.jsonl"
    exit_code = main(
        [
            "view",
            str(indexed_vcf),
            "-r",
            "scaffold_1:4300-5000",
            "--output-format",
            "jsonl",
            "--output-file",
            str(out_file),
            "--hap-prefix",
            "H",
            "--hap-pad",
            "3",
        ]
    )

    assert exit_code == 0
    payload = json.loads(out_file.read_text(encoding="utf-8").strip())
    assert payload["haplotype_label"] == {"prefix": "H", "pad": 3}
    assert payload["haplotypes"][0]["id"] == "H001"
