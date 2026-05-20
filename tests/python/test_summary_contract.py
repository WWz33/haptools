from __future__ import annotations

from haplokit.summary_contract import (
    build_hap_label_map,
    hap_samples,
    hap_states,
    hap_summary_states,
    with_population_breakdown,
)


def test_build_hap_label_map_prefers_contract_id() -> None:
    summary_row = {
        "haplotypes": [
            {"hap": "0/0|1/1", "id": "Hap01"},
            {"hap": "1/1|0/0"},
        ]
    }

    assert build_hap_label_map(summary_row) == {
        "0/0|1/1": "Hap01",
        "1/1|0/0": "Hap02",
    }


def test_hap_states_converts_exact_genotype_states_to_display_alleles() -> None:
    sites = [
        {"allele": "A/T"},
        {"allele": "G/C,A"},
    ]

    assert hap_states("0/0|2/2", sites, "exact") == ["A", "A"]
    assert hap_states("A001", sites, "max-diff") == ["A001"]


def test_hap_summary_states_and_samples_fall_back_to_detail_rows() -> None:
    hap = {"hap": "0/0|1/1"}
    detail_row = {
        "accessions": [
            {"sample": "S1", "hap": "0/0|1/1"},
            {"sample": "S2", "hap": "1/1|0/0"},
        ]
    }
    sites = [{"allele": "A/T"}, {"allele": "G/C"}]

    assert hap_summary_states(hap, sites, "exact") == ["A", "C"]
    assert hap_samples(hap, detail_row) == ["S1"]


def test_summary_contract_uses_backend_states_and_samples_when_present() -> None:
    hap = {
        "hap": "0/0|1/1",
        "states": ["A", "C"],
        "samples": ["S1", "S2"],
    }

    assert hap_summary_states(hap, [], "max-diff") == ["A", "C"]
    assert hap_samples(hap, {"accessions": []}) == ["S1", "S2"]


def test_with_population_breakdown_adds_frequency_labels(tmp_path) -> None:
    pop_file = tmp_path / "pop.tsv"
    pop_file.write_text("S1\tPopA\nS2\tPopA\nS3\tPopB\n", encoding="utf-8")
    haplotypes = [
        {"id": "Hap01", "samples": ["S1", "S3"]},
        {"id": "Hap02", "samples": ["S2"]},
    ]

    enriched = with_population_breakdown(haplotypes, str(pop_file))

    assert enriched[0]["populations"] == [
        {
            "population": "PopA",
            "count": 1,
            "total": 2,
            "frequency": 0.5,
            "frequency_label": "1/2",
        },
        {
            "population": "PopB",
            "count": 1,
            "total": 1,
            "frequency": 1.0,
            "frequency_label": "1/1",
        },
    ]
    assert enriched[1]["populations"][0]["frequency_label"] == "1/2"
