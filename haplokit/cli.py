from __future__ import annotations

import argparse
import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from haplokit.summary_contract import (
    build_hap_label_map,
    hap_samples,
    hap_states,
    hap_summary_states,
    with_population_breakdown,
)


def _region_value(value: str) -> str:
    if ":" not in value:
        raise argparse.ArgumentTypeError("region must look like chr:start-end or chr:pos")
    chrom, coords = value.split(":", 1)
    if not chrom or not coords:
        raise argparse.ArgumentTypeError("region must look like chr:start-end or chr:pos")
    if "-" in coords:
        start, end = coords.split("-", 1)
        if not start.isdigit() or not end.isdigit():
            raise argparse.ArgumentTypeError("region range must be numeric")
    else:
        if not coords.isdigit():
            raise argparse.ArgumentTypeError("site position must be numeric")
    return value


def _max_diff_value(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("max diff must be a float in [0, 1]") from exc
    if parsed < 0 or parsed > 1:
        raise argparse.ArgumentTypeError("max diff must be a float in [0, 1]")
    return parsed


def _hap_pad_value(value: str) -> int:
    if not value.isascii() or not value.isdecimal():
        raise argparse.ArgumentTypeError("hap pad must be a positive integer")
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("hap pad must be a positive integer")
    return parsed


@dataclass(frozen=True)
class Selector:
    payload: dict[str, object]
    region: str


class HaolokitArgumentParser(argparse.ArgumentParser):
    def parse_args(self, args: list[str] | None = None, namespace=None):
        ns = super().parse_args(args=args, namespace=namespace)
        self._validate(ns)
        return ns

    def _validate(self, ns) -> None:
        if getattr(ns, "command", None) != "view":
            return
        has_region = bool(ns.region)
        has_regions_file = bool(ns.regions_file)
        if not has_region and not has_regions_file:
            self.error("one of -r/--region or -R/--regions-file is required")
        if has_region and has_regions_file:
            self.error("-r/--region and -R/--regions-file are mutually exclusive")
        if has_region:
            coords = ns.region.split(":", 1)[1]
            inferred_by = "site" if "-" not in coords else "region"
            if ns.by != "auto" and ns.by != inferred_by:
                self.error(f"--by {ns.by} conflicts with -r selector semantics ({inferred_by})")
            ns.by = inferred_by
            return

        if ns.by == "site":
            self.error("--by site is only valid with -r chr:pos")
        ns.by = "region"


def build_parser() -> HaolokitArgumentParser:
    parser = HaolokitArgumentParser(prog="haplokit")
    subparsers = parser.add_subparsers(dest="command")

    view = subparsers.add_parser("view")
    view.add_argument("input_vcf", nargs="?", default=None)
    view.add_argument("-r", "--region", dest="region", type=_region_value)
    view.add_argument("-R", "--regions-file", dest="regions_file")
    view.add_argument("-S", "--samples-file", dest="samples_file")
    view.add_argument("--by", choices=["auto", "region", "site"], default="auto")
    view.add_argument("--impute", action="store_true")
    view.add_argument("-g", "--gff3", "--gff")
    view.add_argument("--output", dest="output_mode", choices=["summary", "detail"], default="summary")
    view.add_argument("--output-format", choices=["tsv", "jsonl"], default="tsv")
    view.add_argument("--output-file")
    view.add_argument("--plot", action="store_true")
    view.add_argument("--plot-format", choices=["png", "pdf", "svg", "tiff"], default="png")
    view.add_argument("-p", "--population", dest="population_file")
    view.add_argument("--max-diff", type=_max_diff_value)
    view.add_argument("--geo", dest="geo_file")
    view.add_argument("--network", action="store_true")
    view.add_argument("--network-method", choices=["tcs", "msn", "mjn"], default="tcs")
    view.add_argument("--hap-prefix", default="Hap")
    view.add_argument("--hap-pad", type=_hap_pad_value, default=2)

    return parser


def _selector_payload_from_region(region: str) -> tuple[dict[str, object], str]:
    chrom, coords = region.split(":", 1)
    if "-" in coords:
        start, end = coords.split("-", 1)
        return (
            {"type": "region", "chrom": chrom, "start": int(start), "end": int(end)},
            f"{chrom}:{start}-{end}",
        )
    pos = int(coords)
    return ({"type": "site", "chrom": chrom, "pos": pos}, f"{chrom}:{pos}-{pos}")


def _selectors_from_args(args) -> list[Selector]:
    if args.region:
        payload, region = _selector_payload_from_region(args.region)
        return [Selector(payload=payload, region=region)]

    selectors: list[Selector] = []
    for idx, line in enumerate(Path(args.regions_file).read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError("BED rows must have at least 3 columns")
        chrom, start, end = fields[:3]
        record_id = fields[3] if len(fields) > 3 and fields[3] else f"row{idx}"
        selectors.append(
            Selector(
                payload={
                    "type": "bed-record",
                    "record_id": record_id,
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                },
                region=f"{chrom}:{start}-{end}",
            )
        )
    return selectors


def _write_jsonl(rows: Iterable[dict[str, object]], output_file: str | None) -> None:
    data = "\n".join(json.dumps(row, sort_keys=True) for row in rows) + "\n"
    if output_file:
        path = Path(output_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(data, encoding="utf-8")
    else:
        print(data, end="")


def _cpp_backend_path() -> Path:
    env_path = os.environ.get("HAPLOKIT_CPP_BIN")
    candidates: list[Path] = []
    if env_path:
        candidates.append(Path(env_path))

    package_bin = Path(__file__).resolve().parent / "_bin" / "haplokit_cpp"
    candidates.append(package_bin)

    repo_root = Path(__file__).resolve().parents[1]
    candidates.extend(
        [
            repo_root / "build-wsl" / "haplokit_cpp",
            repo_root / "build" / "haplokit_cpp",
        ]
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate

    repo_root = Path(__file__).resolve().parents[1]
    build_dir = repo_root / "build-wsl"
    try:
        subprocess.run(["cmake", "-S", str(repo_root), "-B", str(build_dir)], check=True, cwd=repo_root)
        subprocess.run(["cmake", "--build", str(build_dir), "--clean-first", "-j1"], check=True, cwd=repo_root)
    except Exception:
        pass

    refreshed_candidates = []
    if env_path:
        refreshed_candidates.append(Path(env_path))
    refreshed_candidates.extend([build_dir / "haplokit_cpp", repo_root / "build" / "haplokit_cpp"])
    for candidate in refreshed_candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "haplokit_cpp backend not found; build the C++ target first or set HAPLOKIT_CPP_BIN"
    )


def _run_cpp_view(selector: Selector, args) -> dict[str, object]:
    return _run_cpp_view_mode(selector, args, output_mode=args.output_mode)


def _append_common_args(cmd: list[str], args) -> None:
    if args.samples_file:
        cmd.extend(["--samples-file", str(args.samples_file)])
    if args.impute:
        cmd.append("--impute")
    if args.max_diff is not None:
        cmd.extend(["--max-diff", str(args.max_diff)])
    if args.gff3:
        cmd.extend(["--gff3", str(args.gff3)])
    cmd.extend(["--hap-prefix", str(args.hap_prefix)])
    cmd.extend(["--hap-pad", str(args.hap_pad)])


def _check_backend_result(completed: subprocess.CompletedProcess) -> None:
    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        stdout = completed.stdout.strip()
        detail = stderr or stdout or f"backend exited with code {completed.returncode}"
        raise RuntimeError(detail)


def _run_cpp_view_mode(selector: Selector, args, output_mode: str) -> dict[str, object]:
    cmd = [
        str(_cpp_backend_path()),
        "view-json",
        str(args.input_vcf),
        selector.region,
        "--by",
        args.by,
        "--output",
        output_mode,
    ]
    _append_common_args(cmd, args)

    completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
    _check_backend_result(completed)
    payload = completed.stdout.strip()
    if not payload:
        raise RuntimeError("backend produced no output")
    return json.loads(payload)


def _run_cpp_view_batch(args, output_mode: str) -> list[dict[str, object]]:
    cmd = [
        str(_cpp_backend_path()),
        "view-bed-jsonl",
        str(args.input_vcf),
        str(args.regions_file),
        "--output",
        output_mode,
    ]
    _append_common_args(cmd, args)

    completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
    _check_backend_result(completed)
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    return [json.loads(line) for line in lines]


def _region_slug(selector: Selector) -> str:
    selector_type = selector.payload["type"]
    if selector_type == "site":
        return f'{selector.payload["chrom"]}_{selector.payload["pos"]}_{selector.payload["pos"]}'
    return f'{selector.payload["chrom"]}_{selector.payload["start"]}_{selector.payload["end"]}'


def _sanitized_slug(selector: Selector) -> str:
    return "".join(char if char.isalnum() or char in {"_", "-", "."} else "_" for char in _region_slug(selector))


def _output_base_dir(args) -> Path:
    if not args.output_file:
        return Path.cwd()
    path = Path(args.output_file)
    if path.suffix:
        return path.parent if path.parent != Path("") else Path.cwd()
    return path


def _jsonl_output_path(args) -> str | None:
    if args.output_format != "jsonl":
        return None
    if not args.output_file:
        return None
    path = Path(args.output_file)
    if path.suffix:
        return str(path)
    return str(path / "result.jsonl")


def _output_name_prefix(args) -> str | None:
    if not args.output_file:
        return None
    path = Path(args.output_file)
    if not path.suffix:
        return None
    return path.stem


def _plot_path_for_selector(args, selector: Selector, selector_index: int, selector_count: int) -> Path:
    output_dir = _output_base_dir(args)
    output_dir.mkdir(parents=True, exist_ok=True)
    region_slug = _sanitized_slug(selector)
    fmt = getattr(args, "plot_format", "png")
    path = output_dir / f"{region_slug}.{fmt}"
    if selector_count > 1 and path.exists():
        path = output_dir / f"{region_slug}_{selector_index + 1:03d}.{fmt}"
    return path


def _tsv_paths_for_selector(args, selector: Selector, selector_index: int, selector_count: int) -> tuple[Path, Path]:
    output_dir = _output_base_dir(args)
    output_dir.mkdir(parents=True, exist_ok=True)
    name_prefix = _output_name_prefix(args)
    region_slug = _sanitized_slug(selector)
    if selector_count == 1:
        if name_prefix:
            summary_name = f"{name_prefix}.hap_summary.tsv"
            result_name = f"{name_prefix}.hapresult.tsv"
        else:
            summary_name = "hap_summary.tsv"
            result_name = "hapresult.tsv"
    else:
        if name_prefix:
            summary_name = f"{name_prefix}.hap_summary_{region_slug}.tsv"
            result_name = f"{name_prefix}.hapresult_{region_slug}.tsv"
        else:
            summary_name = f"hap_summary_{region_slug}.tsv"
            result_name = f"hapresult_{region_slug}.tsv"

    summary_path = output_dir / summary_name
    result_path = output_dir / result_name
    if selector_count > 1 and summary_path.exists():
        summary_path = output_dir / f"{summary_path.stem}_{selector_index + 1:03d}{summary_path.suffix}"
    if selector_count > 1 and result_path.exists():
        result_path = output_dir / f"{result_path.stem}_{selector_index + 1:03d}{result_path.suffix}"
    return summary_path, result_path


def _info_cells(site_count: int, annotation: dict[str, object] | None) -> list[str]:
    if site_count <= 0:
        return []
    cells = ["."] * site_count
    if annotation and annotation.get("mode") != "none":
        cells[0] = _annotation_text(annotation)
    return cells


def _write_selector_summary_txt(
    selector: Selector,
    summary_row: dict[str, object],
    detail_row: dict[str, object],
    out_path: Path,
) -> None:
    annotation = summary_row.get("annotation")
    sites = list(summary_row.get("sites", []))
    site_chroms = [str(site["chrom"]) for site in sites]
    site_positions = [str(site["pos"]) for site in sites]
    site_alleles = [str(site["allele"]) for site in sites]
    info_cells = _info_cells(len(sites), annotation if isinstance(annotation, dict) else None)
    lines: list[str] = []
    lines.append("\t".join(["CHR", *site_chroms, "Haplotypes: ", str(summary_row["haplotype_count"])]))
    lines.append("\t".join(["POS", *site_positions, "Individuals: ", str(summary_row["sample_count"])]))
    lines.append("\t".join(["INFO", *info_cells, "Variants: ", str(summary_row["variant_count"])]))
    lines.append("\t".join(["ALLELE", *site_alleles, "Accession", "freq"]))

    for index, item in enumerate(summary_row.get("haplotypes", []), start=1):
        hap_label = str(item.get("id", f"Hap{index:02d}"))
        states = hap_summary_states(item, sites, str(summary_row["grouping_method"]))
        accessions = ";".join(hap_samples(item, detail_row))
        lines.append("\t".join([hap_label, *states, accessions, str(item["count"])]))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_selector_result_txt(
    selector: Selector,
    summary_row: dict[str, object],
    detail_row: dict[str, object],
    out_path: Path,
) -> None:
    annotation = summary_row.get("annotation")
    sites = list(summary_row.get("sites", []))
    site_chroms = [str(site["chrom"]) for site in sites]
    site_positions = [str(site["pos"]) for site in sites]
    site_alleles = [str(site["allele"]) for site in sites]
    info_cells = _info_cells(len(sites), annotation if isinstance(annotation, dict) else None)
    lines: list[str] = []
    lines.append("\t".join(["CHR", *site_chroms, "Haplotypes: ", str(summary_row["haplotype_count"])]))
    lines.append("\t".join(["POS", *site_positions, "Individuals: ", str(summary_row["sample_count"])]))
    lines.append("\t".join(["INFO", *info_cells, "Variants: ", str(summary_row["variant_count"])]))
    lines.append("\t".join(["ALLELE", *site_alleles, "Accession"]))

    hap_label_map = build_hap_label_map(summary_row)
    for item in detail_row.get("accessions", []):
        hap_label = hap_label_map.get(item["hap"], item["hap"])
        summary_item = next(
            (hap for hap in summary_row.get("haplotypes", []) if hap.get("hap") == item["hap"]),
            None,
        )
        states = (
            hap_summary_states(summary_item, sites, str(summary_row["grouping_method"]))
            if isinstance(summary_item, dict)
            else hap_states(item["hap"], sites, str(summary_row["grouping_method"]))
        )
        lines.append("\t".join([hap_label, *states, str(item["sample"])]))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _selector_span(selector: Selector) -> tuple[str, int, int]:
    selector_type = selector.payload["type"]
    chrom = str(selector.payload["chrom"])
    if selector_type == "site":
        pos = int(selector.payload["pos"])
        return chrom, pos, pos
    return chrom, int(selector.payload["start"]), int(selector.payload["end"])


def _annotation_text(annotation: dict[str, object]) -> str:
    mode = annotation.get("mode")
    if mode == "overlap":
        return str(annotation.get("id", "NA"))
    if mode == "nearest":
        return f'{annotation.get("id", "NA")}:{annotation.get("distance", "NA")}'
    return "NA"


def _write_gff_annotation_table(rows: list[dict[str, object]], selectors: list[Selector], args) -> Path:
    output_dir = _output_base_dir(args)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "gff_ann_summary.tsv"
    lines = ["chr\tstart\tend\tann"]
    for selector, row in zip(selectors, rows):
        chrom, start, end = _selector_span(selector)
        annotation = row.get("annotation", {"mode": "none"})
        lines.append(f"{chrom}\t{start}\t{end}\t{_annotation_text(annotation)}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def _compose_row(
    selector: Selector,
    args,
    backend_row: dict[str, object],
) -> dict[str, object]:
    row: dict[str, object] = {
        "grouping_mode": backend_row["grouping_mode"],
        "grouping_method": backend_row["grouping_method"],
        "haplotype_label": backend_row.get(
            "haplotype_label",
            {"prefix": args.hap_prefix, "pad": args.hap_pad},
        ),
        "output_mode": backend_row["output_mode"],
        "imputed_ref": backend_row["imputed_ref"],
        "plot_requested": bool(args.plot),
        "selector": selector.payload,
        "region_label": selector.region,
        "gff3_enabled": bool(args.gff3),
        "max_diff": backend_row["max_diff"],
        "variant_count": backend_row["variant_count"],
        "sample_count": backend_row["sample_count"],
        "haplotype_count": backend_row["haplotype_count"],
        "sites": backend_row.get("sites", []),
    }
    if args.plot:
        row["plot_backend"] = "python"

    if args.output_mode == "summary":
        row["haplotypes"] = with_population_breakdown(list(backend_row["haplotypes"]), args.population_file)
    elif args.output_mode == "detail":
        row["accessions"] = backend_row["accessions"]
    else:  # both
        row["haplotypes"] = with_population_breakdown(list(backend_row.get("haplotypes", [])), args.population_file)
        row["accessions"] = backend_row.get("accessions", [])

    if "annotation" in backend_row:
        row["annotation"] = backend_row["annotation"]

    return row


def _write_plot_artifacts(
    args,
    selectors: list[Selector],
    summary_rows: list[dict[str, object]],
    detail_rows: list[dict[str, object]],
) -> list[dict[str, str]]:
    written: list[dict[str, str]] = []
    for idx, selector in enumerate(selectors):
        summary_row = summary_rows[idx]
        detail_row = detail_rows[idx]

        summary_path, result_path = _tsv_paths_for_selector(args, selector, idx, len(selectors))
        _write_selector_summary_txt(selector, summary_row, detail_row, summary_path)
        _write_selector_result_txt(selector, summary_row, detail_row, result_path)
        artifact_paths: dict[str, str] = {
            "hap_summary_file": str(summary_path.resolve()),
            "hapresult_file": str(result_path.resolve()),
        }
        artifact_paths["summary_file"] = artifact_paths["hap_summary_file"]

        if args.plot:
            from haplokit.plot import plot_hap_table, read_hap_summary_tsv, read_popgroup

            pdf_path = _plot_path_for_selector(args, selector, idx, len(selectors))
            summary_table = read_hap_summary_tsv(summary_path)
            gene_name = ""
            ann = summary_rows[idx].get("annotation")
            if isinstance(ann, dict) and ann.get("mode") != "none":
                gene_name = str(ann.get("id", ""))
            pop_data = None
            if args.population_file:
                pop_data = read_popgroup(args.population_file)
            gff_path = str(args.gff3) if args.gff3 else None
            artifact_paths["plot_file"] = str(plot_hap_table(summary_table, pdf_path, pop_data=pop_data, gff_path=gff_path, title=gene_name, fmt=args.plot_format))

        # Geographic distribution map
        if args.geo_file:
            from haplokit.plot import plot_hap_distribution, read_popgroup

            detail_row = detail_rows[idx]
            summary_row_data = summary_rows[idx]
            haplotypes = list(summary_row_data.get("haplotypes", []))
            hap_names = [str(h.get("id", f"Hap{i:02d}")) for i, h in enumerate(haplotypes, 1)]
            hap_colors = None
            if len(hap_names) > 0:
                from haplokit._palette import allele_palette
                hap_colors = list(allele_palette(hap_names).values())

            # Build sample→haplotype mapping from detail accessions
            sample_hap: dict[str, str] = {}
            for hi, hap in enumerate(haplotypes):
                label = hap_names[hi] if hi < len(hap_names) else str(hap.get("hap", ""))
                for sample in hap_samples(hap, detail_row):
                    sample_hap[sample] = label

            # Read geo file
            geo_samples: list[dict] = []
            with Path(args.geo_file).open("r", encoding="utf-8") as gf:
                import csv
                reader = csv.DictReader(gf, delimiter="\t")
                for row in reader:
                    sid = row.get("ID", "").strip()
                    if sid in sample_hap:
                        geo_samples.append({
                            "lon": float(row.get("longitude", 0)),
                            "lat": float(row.get("latitude", 0)),
                            "hap": sample_hap[sid],
                        })

            if geo_samples:
                geo_path = _plot_path_for_selector(args, selector, idx, len(selectors))
                geo_path = geo_path.with_name(geo_path.stem + "_map" + geo_path.suffix)
                artifact_paths["geo_file"] = str(plot_hap_distribution(
                    geo_samples, hap_names, geo_path,
                    hap_colors=hap_colors, title=gene_name,
                    fmt=args.plot_format,
                ))

        # Haplotype network
        if args.network:
            from haplokit.plot import plot_hap_network, read_popgroup

            summary_row_data = summary_rows[idx]
            detail_row = detail_rows[idx]
            haplotypes = list(summary_row_data.get("haplotypes", []))
            hap_names = [str(h.get("id", f"Hap{i:02d}")) for i, h in enumerate(haplotypes, 1)]
            hap_strings = [str(h.get("pattern", h["hap"])) for h in haplotypes]
            hap_counts = [int(h["count"]) for h in haplotypes]

            # Build population composition per haplotype
            pop_data_net = None
            if args.population_file:
                pop_map = read_popgroup(args.population_file)
                pop_names = sorted(set(pop_map.values()))
                pop_data_net = {}
                for hi, hap in enumerate(haplotypes):
                    label = hap_names[hi] if hi < len(hap_names) else str(hap.get("hap", ""))
                    counts_by_pop: dict[str, int] = {}
                    for s in hap_samples(hap, detail_row):
                        p = pop_map.get(s, "Unknown")
                        counts_by_pop[p] = counts_by_pop.get(p, 0) + 1
                    pop_data_net[label] = [(p, counts_by_pop[p]) for p in pop_names if p in counts_by_pop]

            net_path = _plot_path_for_selector(args, selector, idx, len(selectors))
            net_path = net_path.with_name(net_path.stem + "_network" + net_path.suffix)
            artifact_paths["network_file"] = str(plot_hap_network(
                hap_names, hap_strings, hap_counts, net_path,
                pop_data=pop_data_net, title=gene_name,
                fmt=args.plot_format,
                algorithm=args.network_method,
            ))
        written.append(artifact_paths)

    return written


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command != "view":
        parser.error("a subcommand is required")
    selectors = _selectors_from_args(args)
    if args.output_format == "tsv":
        if args.regions_file:
            summary_rows = _run_cpp_view_batch(args, "summary")
            detail_rows = _run_cpp_view_batch(args, "detail")
        else:
            both_rows = [_run_cpp_view_mode(selector, args, "both") for selector in selectors]
            summary_rows = both_rows
            detail_rows = both_rows

        if len(summary_rows) != len(selectors) or len(detail_rows) != len(selectors):
            raise RuntimeError("backend row count did not match selector count")

        _write_plot_artifacts(args, selectors, summary_rows, detail_rows)

        if args.gff3:
            rows = [{"annotation": summary_rows[idx].get("annotation", {"mode": "none"})} for idx in range(len(selectors))]
            _write_gff_annotation_table(rows, selectors, args)
        return 0

    if args.regions_file:
        backend_rows = _run_cpp_view_batch(args, args.output_mode)
        if len(backend_rows) != len(selectors):
            raise RuntimeError("backend row count did not match selector count")
    else:
        backend_rows = [_run_cpp_view(selector, args) for selector in selectors]

    rows = [
        _compose_row(
            selector,
            args,
            backend_row=backend_rows[index],
        )
        for index, selector in enumerate(selectors)
    ]

    if args.plot:
        if args.regions_file:
            both_rows = _run_cpp_view_batch(args, "both")
            summary_rows = both_rows
            detail_rows = both_rows
        else:
            both_rows = [_run_cpp_view_mode(selector, args, "both") for selector in selectors]
            summary_rows = both_rows
            detail_rows = both_rows

        written = _write_plot_artifacts(args, selectors, summary_rows, detail_rows)
        if len(written) != len(rows):
            raise RuntimeError("plot artifact count did not match selector count")
        for idx, row in enumerate(rows):
            row.update(written[idx])

    if args.gff3:
        gff_summary = _write_gff_annotation_table(rows, selectors, args)
        for row in rows:
            row["gff_summary_file"] = str(gff_summary.resolve())
    _write_jsonl(rows, _jsonl_output_path(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
