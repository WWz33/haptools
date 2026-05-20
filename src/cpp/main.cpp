#include <iostream>
#include <optional>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>

#include "gff_annotator.h"
#include "selector.h"
#include "vcf_reader.h"
#include "view_backend.h"

namespace {

struct ParsedViewJsonCommand {
    std::string input_vcf;
    haplokit::Region region;
    std::vector<std::string> samples;
    haplokit::ViewOptions options;
    std::string gff3_path;
};

double parse_max_diff(const std::string& value) {
    const double parsed = std::stod(value);
    if (parsed < 0.0 || parsed > 1.0) {
        throw std::runtime_error("max diff must be in [0,1]");
    }
    return parsed;
}

int parse_hap_pad(const std::string& value) {
    const int parsed = std::stoi(value);
    if (parsed < 1) {
        throw std::runtime_error("hap pad must be a positive integer");
    }
    return parsed;
}

ParsedViewJsonCommand parse_view_json_command(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp view-json <vcf> <region> [--by region|site] [--samples-file path] "
            "[--impute] [--output summary|detail] [--max-diff x]");
    }

    ParsedViewJsonCommand parsed;
    parsed.input_vcf = argv[2];
    parsed.region = haplokit::parse_region(argv[3]);

    for (int idx = 4; idx < argc; ++idx) {
        const std::string arg = argv[idx];
        if (arg == "--by") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--by requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "region") {
                parsed.options.by = haplokit::GroupBy::Region;
            } else if (value == "site") {
                parsed.options.by = haplokit::GroupBy::Site;
            } else {
                throw std::runtime_error("unsupported --by value: " + value);
            }
            continue;
        }
        if (arg == "--samples-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--samples-file requires a path");
            }
            parsed.samples = haplokit::load_sample_list(argv[++idx]);
            continue;
        }
        if (arg == "--impute") {
            parsed.options.impute = true;
            continue;
        }
        if (arg == "--output") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--output requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "summary") {
                parsed.options.output_mode = haplokit::OutputMode::Summary;
            } else if (value == "detail") {
                parsed.options.output_mode = haplokit::OutputMode::Detail;
            } else if (value == "both") {
                parsed.options.output_mode = haplokit::OutputMode::Both;
            } else {
                throw std::runtime_error("unsupported --output value: " + value);
            }
            continue;
        }
        if (arg == "--max-diff") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--max-diff requires a value");
            }
            parsed.options.max_diff = parse_max_diff(argv[++idx]);
            continue;
        }
        if (arg == "--gff3") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--gff3 requires a path");
            }
            parsed.gff3_path = argv[++idx];
            continue;
        }
        if (arg == "--hap-prefix") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-prefix requires a value");
            }
            parsed.options.hap_prefix = argv[++idx];
            continue;
        }
        if (arg == "--hap-pad") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-pad requires a value");
            }
            parsed.options.hap_pad = parse_hap_pad(argv[++idx]);
            continue;
        }
        throw std::runtime_error("unsupported argument: " + arg);
    }

    if (parsed.options.by == haplokit::GroupBy::Site && !haplokit::is_site_region(parsed.region)) {
        throw std::runtime_error("--by site requires a single-position region");
    }

    return parsed;
}

haplokit::ViewOptions parse_common_view_options(
    int argc,
    char** argv,
    int start_index,
    std::vector<std::string>* samples,
    std::string* gff3_path) {
    haplokit::ViewOptions options;
    if (samples == nullptr) {
        throw std::runtime_error("samples output is required");
    }
    for (int idx = start_index; idx < argc; ++idx) {
        const std::string arg = argv[idx];
        if (arg == "--samples-file") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--samples-file requires a path");
            }
            *samples = haplokit::load_sample_list(argv[++idx]);
            continue;
        }
        if (arg == "--impute") {
            options.impute = true;
            continue;
        }
        if (arg == "--output") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--output requires a value");
            }
            const std::string value = argv[++idx];
            if (value == "summary") {
                options.output_mode = haplokit::OutputMode::Summary;
            } else if (value == "detail") {
                options.output_mode = haplokit::OutputMode::Detail;
            } else {
                throw std::runtime_error("unsupported --output value: " + value);
            }
            continue;
        }
        if (arg == "--max-diff") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--max-diff requires a value");
            }
            options.max_diff = parse_max_diff(argv[++idx]);
            continue;
        }
        if (arg == "--gff3") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--gff3 requires a path");
            }
            if (gff3_path) {
                *gff3_path = argv[++idx];
            }
            continue;
        }
        if (arg == "--hap-prefix") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-prefix requires a value");
            }
            options.hap_prefix = argv[++idx];
            continue;
        }
        if (arg == "--hap-pad") {
            if (idx + 1 >= argc) {
                throw std::runtime_error("--hap-pad requires a value");
            }
            options.hap_pad = parse_hap_pad(argv[++idx]);
            continue;
        }
        throw std::runtime_error("unsupported argument: " + arg);
    }
    return options;
}

std::vector<haplokit::Region> load_bed_regions(const std::string& path) {
    std::ifstream handle(path);
    if (!handle) {
        throw std::runtime_error("failed to open regions file: " + path);
    }

    std::vector<haplokit::Region> regions;
    std::string line;
    while (std::getline(handle, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const auto end = line.find('\t', start);
            if (end == std::string::npos) {
                fields.push_back(line.substr(start));
                break;
            }
            fields.push_back(line.substr(start, end - start));
            start = end + 1;
        }
        if (fields.size() < 3) {
            throw std::runtime_error("BED rows must have at least 3 columns");
        }

        haplokit::Region region;
        region.chrom = fields[0];
        region.start = std::stoi(fields[1]);
        region.end = std::stoi(fields[2]);
        regions.push_back(region);
    }
    return regions;
}

int run_view_json(int argc, char** argv) {
    const auto parsed = parse_view_json_command(argc, argv);

    std::optional<haplokit::GffAnnotator> annotator;
    if (!parsed.gff3_path.empty()) {
        haplokit::GffAnnotator ann;
        if (!ann.load(parsed.gff3_path)) {
            throw std::runtime_error("failed to parse GFF3: " + parsed.gff3_path);
        }
        annotator = std::move(ann);
    }

    haplokit::VcfReader reader(parsed.input_vcf);
    const auto data = reader.fetch(parsed.region, parsed.samples);
    auto result = haplokit::build_view_result(data, parsed.options);

    if (annotator.has_value()) {
        result.annotation = annotator->annotate(parsed.region.chrom, parsed.region.start, parsed.region.end);
    }

    std::cout << haplokit::serialize_view_result_json(result) << "\n";
    return 0;
}

int run_view_bed_jsonl(int argc, char** argv) {
    if (argc < 4) {
        throw std::runtime_error(
            "usage: haplokit_cpp view-bed-jsonl <vcf> <regions.bed> [--samples-file path] [--impute] "
            "[--output summary|detail] [--max-diff x] [--gff3 path]");
    }

    const std::string input_vcf = argv[2];
    const std::string regions_file = argv[3];
    std::vector<std::string> samples;
    std::string gff3_path;
    auto options = parse_common_view_options(argc, argv, 4, &samples, &gff3_path);
    options.by = haplokit::GroupBy::Region;

    std::optional<haplokit::GffAnnotator> annotator;
    if (!gff3_path.empty()) {
        haplokit::GffAnnotator ann;
        if (!ann.load(gff3_path)) {
            throw std::runtime_error("failed to parse GFF3: " + gff3_path);
        }
        annotator = std::move(ann);
    }

    haplokit::VcfReader reader(input_vcf);
    const auto regions = load_bed_regions(regions_file);
    for (std::size_t idx = 0; idx < regions.size(); ++idx) {
        const auto data = reader.fetch(regions[idx], samples);
        auto result = haplokit::build_view_result(data, options);

        if (annotator.has_value()) {
            // BED is 0-based half-open; GFF annotation needs 1-based inclusive
            int64_t gff_start = regions[idx].start + 1;
            int64_t gff_end = regions[idx].end;
            result.annotation = annotator->annotate(regions[idx].chrom, gff_start, gff_end);
        }

        std::cout << haplokit::serialize_view_result_json(result) << "\n";
    }
    return 0;
}

int run_debug_fetch(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: haplokit_cpp <vcf> <region>\n";
        return 1;
    }

    const haplokit::Region region = haplokit::parse_region(argv[2]);
    haplokit::VcfReader reader(argv[1]);
    const auto data = reader.fetch(region);
    std::cout << "opened " << argv[1] << " for " << region.chrom << ":" << region.start << "-" << region.end << "\n";
    std::cout << "samples " << data.samples.size() << "\n";
    std::cout << "variants " << data.variants.size() << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 2 && std::string(argv[1]) == "view-json") {
            return run_view_json(argc, argv);
        }
        if (argc >= 2 && std::string(argv[1]) == "view-bed-jsonl") {
            return run_view_bed_jsonl(argc, argv);
        }
        return run_debug_fetch(argc, argv);
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 2;
    }
}
