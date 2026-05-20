#pragma once

#include <optional>
#include <string>
#include <vector>

#include "gff_annotator.h"
#include "vcf_reader.h"

namespace haplokit {

enum class GroupBy {
    Region,
    Site,
};

enum class OutputMode {
    Summary,
    Detail,
    Both,
};

struct ViewOptions {
    GroupBy by = GroupBy::Region;
    bool impute = false;
    std::optional<double> max_diff;
    OutputMode output_mode = OutputMode::Summary;
    std::string hap_prefix = "Hap";
    int hap_pad = 2;
};

struct HaplotypeSummaryRow {
    std::string id;
    std::string hap;
    std::vector<std::string> states;
    std::vector<std::string> samples;
    int count = 0;
    int total = 0;
    double frequency = 0.0;
    std::string frequency_label;
};

struct HaplotypeDetailRow {
    std::string sample;
    std::string hap;
};

struct SiteRow {
    std::string chrom;
    int pos = 0;
    std::string allele;
};

struct ViewResult {
    std::string grouping_mode;
    std::string grouping_method;
    std::string output_mode;
    bool imputed_ref = false;
    std::optional<double> max_diff;
    std::string hap_prefix = "Hap";
    int hap_pad = 2;
    int variant_count = 0;
    int sample_count = 0;
    int haplotype_count = 0;
    std::vector<SiteRow> sites;
    std::vector<HaplotypeSummaryRow> haplotypes;
    std::vector<HaplotypeDetailRow> accessions;
    std::optional<GeneAnnotation> annotation;
};

ViewResult build_view_result(const RegionData& data, const ViewOptions& options);
std::string serialize_view_result_json(const ViewResult& result);
std::vector<std::string> load_sample_list(const std::string& path);

}  // namespace haplokit
