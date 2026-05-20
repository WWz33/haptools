#include <iostream>
#include <stdexcept>

#include "gff_annotator.h"
#include "selector.h"
#include "vcf_reader.h"
#include "view_backend.h"

namespace {

int run() {
    haplokit::GffAnnotator annotator;
    if (!annotator.load("data/annotation.gff")) {
        std::cerr << "failed to load GFF fixture\n";
        return 1;
    }
    const auto gene = annotator.find_gene("test1G0387");
    if (!gene.has_value() || gene->mode != "gene" || gene->seqid != "scaffold_1" ||
        gene->start != 4300 || gene->end != 7910 || gene->strand != '+') {
        std::cerr << "expected gffsub-backed gene lookup for test1G0387\n";
        return 1;
    }
    if (annotator.find_gene("missing_gene").has_value()) {
        std::cerr << "unexpected gene lookup for missing_gene\n";
        return 1;
    }

    haplokit::VcfReader reader("data/var.sorted.vcf.gz");

    const auto region_data = reader.fetch(haplokit::parse_region("scaffold_1:4300-5000"));
    haplokit::ViewOptions strict_region;
    strict_region.by = haplokit::GroupBy::Region;
    strict_region.output_mode = haplokit::OutputMode::Summary;
    const auto summary = haplokit::build_view_result(region_data, strict_region);
    if (summary.grouping_mode != "strict-region") {
        std::cerr << "unexpected grouping mode: " << summary.grouping_mode << "\n";
        return 1;
    }
    if (summary.variant_count != 5 || summary.sample_count != 37) {
        std::cerr << "unexpected counts: variants=" << summary.variant_count
                  << " samples=" << summary.sample_count << "\n";
        return 1;
    }
    if (summary.haplotype_count <= 0 || summary.haplotypes.empty()) {
        std::cerr << "expected non-empty summary haplotypes\n";
        return 1;
    }
    if (summary.sites.size() != 5U || summary.sites.front().pos != 4300) {
        std::cerr << "expected site metadata for summary plot\n";
        return 1;
    }
    const auto& first_hap = summary.haplotypes.front();
    if (summary.hap_prefix != "Hap" || summary.hap_pad != 2 || first_hap.id != "Hap01") {
        std::cerr << "unexpected default haplotype label contract\n";
        return 1;
    }
    if (first_hap.hap.empty() || first_hap.states.empty() || first_hap.samples.empty()) {
        std::cerr << "expected summary haplotype to include pattern, states, and samples\n";
        return 1;
    }
    if (first_hap.total != summary.sample_count || first_hap.frequency <= 0.0 ||
        first_hap.frequency_label != std::to_string(first_hap.count) + "/" + std::to_string(summary.sample_count)) {
        std::cerr << "unexpected summary haplotype count/frequency contract\n";
        return 1;
    }

    haplokit::ViewOptions custom_labels = strict_region;
    custom_labels.hap_prefix = "H";
    custom_labels.hap_pad = 3;
    const auto custom_summary = haplokit::build_view_result(region_data, custom_labels);
    if (custom_summary.haplotype_count <= 0 || custom_summary.haplotypes.front().id != "H001") {
        std::cerr << "custom haplotype labels were not applied\n";
        return 1;
    }

    haplokit::ViewOptions strict_site;
    strict_site.by = haplokit::GroupBy::Site;
    strict_site.output_mode = haplokit::OutputMode::Detail;
    const auto site_data = reader.fetch(haplokit::parse_region("scaffold_1:4300"));
    const auto detail = haplokit::build_view_result(site_data, strict_site);
    if (detail.grouping_mode != "strict-site") {
        std::cerr << "unexpected site grouping mode: " << detail.grouping_mode << "\n";
        return 1;
    }
    if (detail.variant_count != 1 || detail.accessions.empty()) {
        std::cerr << "expected one site variant and non-empty detail rows\n";
        return 1;
    }

    haplokit::ViewOptions approx_region;
    approx_region.by = haplokit::GroupBy::Region;
    approx_region.output_mode = haplokit::OutputMode::Summary;
    approx_region.max_diff = 0.2;
    const auto approx = haplokit::build_view_result(region_data, approx_region);
    if (approx.grouping_mode != "approx-region" || approx.grouping_method != "max-diff") {
        std::cerr << "unexpected approx metadata\n";
        return 1;
    }

    const auto payload = haplokit::serialize_view_result_json(summary);
    if (payload.find("\"grouping_mode\":\"strict-region\"") == std::string::npos) {
        std::cerr << "serialized payload missing grouping_mode\n";
        return 1;
    }
    if (payload.find("\"haplotypes\"") == std::string::npos) {
        std::cerr << "serialized payload missing haplotypes\n";
        return 1;
    }
    if (payload.find("\"haplotype_label\":{\"pad\":2,\"prefix\":\"Hap\"}") == std::string::npos ||
        payload.find("\"frequency_label\"") == std::string::npos ||
        payload.find("\"id\":\"Hap01\"") == std::string::npos ||
        payload.find("\"pattern\"") == std::string::npos ||
        payload.find("\"samples\"") == std::string::npos ||
        payload.find("\"states\"") == std::string::npos) {
        std::cerr << "serialized payload missing summary contract fields\n";
        return 1;
    }
    if (payload.find("\"sites\"") == std::string::npos) {
        std::cerr << "serialized payload missing sites\n";
        return 1;
    }

    return 0;
}

}  // namespace

int main() {
    try {
        return run();
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 2;
    }
}
