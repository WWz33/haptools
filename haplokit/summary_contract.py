from __future__ import annotations

from pathlib import Path

__all__ = [
    "build_hap_label_map",
    "hap_samples",
    "hap_states",
    "hap_summary_states",
    "read_popgroup_file",
    "with_population_breakdown",
]


def _state_to_label(state: str, allele: str) -> str:
    if "/" not in state or "/" not in allele:
        return state
    ref, alt = allele.split("/", 1)
    alts = [item for item in alt.split(",") if item]
    left, _right = state.split("/", 1)
    allele_index = int(left)
    if allele_index == 0:
        return ref
    alt_index = allele_index - 1
    if 0 <= alt_index < len(alts):
        return alts[alt_index]
    return state


def build_hap_label_map(summary_row: dict[str, object]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for index, item in enumerate(summary_row.get("haplotypes", []), start=1):
        if isinstance(item, dict):
            mapping[str(item["hap"])] = str(item.get("id", f"Hap{index:02d}"))
    return mapping


def hap_states(hap_value: str, sites: list[dict[str, object]], grouping_method: str) -> list[str]:
    if grouping_method == "exact" and sites:
        states = hap_value.split("|")
        labels: list[str] = []
        for idx, state in enumerate(states):
            allele = str(sites[idx]["allele"]) if idx < len(sites) else ""
            labels.append(_state_to_label(state, allele) if allele else state)
        return labels
    return [hap_value]


def hap_summary_states(item: dict[str, object], sites: list[dict[str, object]], grouping_method: str) -> list[str]:
    states = item.get("states")
    if isinstance(states, list):
        return [str(state) for state in states]
    return hap_states(str(item["hap"]), sites, grouping_method)


def hap_samples(item: dict[str, object], detail_row: dict[str, object]) -> list[str]:
    samples = item.get("samples")
    if isinstance(samples, list):
        return [str(sample) for sample in samples]
    return [
        str(accession["sample"])
        for accession in detail_row.get("accessions", [])
        if isinstance(accession, dict) and accession.get("hap") == item.get("hap")
    ]


def read_popgroup_file(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    pop: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) >= 2 and fields[0].strip():
                pop[fields[0].strip()] = fields[1].strip()
    return pop


def with_population_breakdown(
    haplotypes: list[dict[str, object]],
    population_file: str | None,
) -> list[dict[str, object]]:
    pop_map = read_popgroup_file(population_file)
    if not pop_map:
        return haplotypes
    population_totals: dict[str, int] = {}
    for population in pop_map.values():
        population_totals[population] = population_totals.get(population, 0) + 1

    enriched: list[dict[str, object]] = []
    for hap in haplotypes:
        item = dict(hap)
        counts: dict[str, int] = {}
        for sample in item.get("samples", []):
            population = pop_map.get(str(sample))
            if population:
                counts[population] = counts.get(population, 0) + 1
        item["populations"] = [
            {
                "population": population,
                "count": count,
                "total": population_totals[population],
                "frequency": count / population_totals[population] if population_totals[population] else 0,
                "frequency_label": f"{count}/{population_totals[population]}",
            }
            for population, count in sorted(counts.items())
        ]
        enriched.append(item)
    return enriched
