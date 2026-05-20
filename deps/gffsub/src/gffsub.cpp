#include "gff3.hpp"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <getopt.h>

using namespace gffsub;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.gff3> [options]\n"
        << "\n"
        << "Input/Region Options:\n"
        << "  -r, --region CHR:START-END\n"
        << "      Extract features overlapping the specified genomic region.\n"
        << "      Coordinates are 1-based and inclusive (GFF format).\n"
        << "      Example: -r chr1:1000000-2000000\n"
        << "\n"
        << "  -b, --bed FILE\n"
        << "      Extract features using genomic regions from a BED file.\n"
        << "      BED files use 0-based half-open coordinates, automatically\n"
        << "      converted to 1-based for internal processing.\n"
        << "\n"
        << "Feature Filter Options:\n"
        << "  -f, --feature TYPE\n"
        << "      Filter features by type (3rd column in GFF/GTF).\n"
        << "      Examples: gene, mRNA, exon, CDS, transcript\n"
        << "\n"
        << "  -L, --longest\n"
        << "      Keep only the longest transcript isoform for each gene.\n"
        << "      Per-gene comparison (AGAT logic): if gene has CDS isoforms,\n"
        << "      only compare by CDS length; otherwise compare by exon length.\n"
        << "\n"
        << "  -@, --threads N\n"
        << "      Number of threads for parallel processing (default: 1).\n"
        << "      Currently used with --longest for multi-chromosome parallelization.\n"
        << "\n"
        << "Output Options:\n"
        << "  -t, --output-format FMT\n"
        << "      Output format. Choices: gff3, gtf2, gtf3, bed\n"
        << "      gff3  - GFF3 format (##gff-version 3)\n"
        << "      gtf2  - GTF2 format (##gtf-version 2)\n"
        << "      gtf3  - GTF3/Ensembl format (##gtf-version 2.2.1)\n"
        << "      bed   - BED format (0-based half-open coordinates)\n"
        << "      Default: gff3\n"
        << "\n"
        << "  -o, --output FILE\n"
        << "      Output file path. If not specified, writes to stdout.\n"
        << "\n"
        << "  -h, --help\n"
        << "      Display this help message.\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -f gene\n"
        << "  " << prog << " annotation.gff3 --bed regions.bed -f exon\n"
        << "  " << prog << " annotation.gff3 --longest\n"
        << "  " << prog << " annotation.gff3 --longest -@ 6\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -t gtf3 -o out.gtf\n";
}

static void query_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " query <input.gff3> [options]\n"
        << "\n"
        << "Query Options:\n"
        << "  --id ID                 Query a feature by ID.\n"
        << "  --name NAME             Query a gene by Name/Alias/gene_id/locus_tag/ID.\n"
        << "  --id-list FILE          Query one feature ID per line.\n"
        << "  --region CHR:START-END  Query features overlapping a 1-based inclusive region.\n"
        << "  --type TYPE             Restrict query output by feature type.\n"
        << "  --attr KEY=VALUE        Query features by an exact GFF3 attribute value.\n"
        << "  --include-children      Include descendants of matched IDs.\n"
        << "  --summary-format FMT    Output query summary instead of GFF3. Choices: tsv, json.\n"
        << "  -h, --help              Display this help message.\n";
}

static void window_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " window <input.gff3> [options]\n"
        << "\n"
        << "Window Options:\n"
        << "  --id ID                 Target feature ID or gene lookup key.\n"
        << "  --upstream N            Bases to add upstream of the target (default: 0).\n"
        << "  --downstream N          Bases to add downstream of the target (default: 0).\n"
        << "  --strand-aware          Interpret upstream/downstream by feature strand.\n"
        << "  -h, --help              Display this help message.\n";
}

static void qc_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " qc <input.gff3>\n"
        << "\n"
        << "QC checks:\n"
        << "  duplicate_id      Repeated ID attributes.\n"
        << "  invalid_range     start greater than end.\n"
        << "  missing_parent    Parent points to an absent ID.\n"
        << "  child_outside_parent  Child coordinates outside parent coordinates.\n";
}

struct SummaryRow {
    std::string query_id;
    std::string matched_id;
    std::string matched_by;
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::string type;
    std::string parent_id;
    size_t child_count = 0;
    size_t transcript_count = 0;
    size_t exon_count = 0;
    int64_t cds_length = 0;
    std::string status;
};

static bool append_unique(GffData& out, std::unordered_set<int>& seen, const GffRecord& rec) {
    if (!seen.insert(rec.line_idx).second) {
        return false;
    }
    out.append(rec);
    return true;
}

static std::string record_id(const GffRecord& rec) {
    if (rec.id) return *rec.id;
    if (rec.gene_id) return *rec.gene_id;
    if (rec.transcript_id) return *rec.transcript_id;
    return "";
}

static void add_feature_counts(SummaryRow& row, const std::vector<GffRecord>& records) {
    for (const auto& rec : records) {
        if (rec.type == "mRNA" || rec.type == "transcript") {
            ++row.transcript_count;
        } else if (rec.type == "exon") {
            ++row.exon_count;
        } else if (rec.type == "CDS") {
            row.cds_length += rec.end - rec.start + 1;
        }
    }
}

static SummaryRow make_summary_row(const AnnotationIndex& index,
                                   const std::string& query_id,
                                   const std::string& matched_by,
                                   const GffRecord& rec) {
    SummaryRow row;
    row.query_id = query_id;
    row.matched_id = record_id(rec);
    row.matched_by = matched_by;
    row.seqid = rec.seqid;
    row.start = rec.start;
    row.end = rec.end;
    row.strand = rec.strand;
    row.type = rec.type;
    row.parent_id = rec.parent_id.value_or("");
    row.status = "found";

    if (rec.id) {
        const auto children = index.children_of(*rec.id);
        row.child_count = children.size();
        const auto model = index.gene_model(*rec.id);
        if (model) {
            add_feature_counts(row, model->records);
        } else {
            add_feature_counts(row, children);
        }
    }

    return row;
}

static SummaryRow make_not_found_row(const std::string& query_id, const std::string& matched_by) {
    SummaryRow row;
    row.query_id = query_id;
    row.matched_by = matched_by;
    row.status = "not_found";
    return row;
}

static bool contains_record(const std::vector<GffRecord>& records, int line_idx) {
    for (const auto& rec : records) {
        if (rec.line_idx == line_idx) {
            return true;
        }
    }
    return false;
}

static std::string infer_gene_match_key(const AnnotationIndex& index, const std::string& query, const GffRecord& rec) {
    if (rec.id && *rec.id == query) {
        return "ID";
    }
    if (rec.gene_id && *rec.gene_id == query) {
        return "gene_id";
    }
    for (const char* key : {"Name", "locus_tag", "Alias", "Dbxref"}) {
        if (contains_record(index.with_attribute(key, query), rec.line_idx)) {
            return key;
        }
    }
    return "name";
}

static std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        switch (ch) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default: out << ch; break;
        }
    }
    return out.str();
}

static std::vector<std::string> attribute_values(const std::string& attrs, const std::string& key) {
    std::vector<std::string> values;
    size_t pos = 0;
    while (pos < attrs.size()) {
        const size_t key_end = attrs.find('=', pos);
        if (key_end == std::string::npos) {
            break;
        }

        const std::string found_key = attrs.substr(pos, key_end - pos);
        size_t value_end = attrs.find(';', key_end + 1);
        if (value_end == std::string::npos) {
            value_end = attrs.size();
        }

        if (found_key == key) {
            size_t value_start = key_end + 1;
            while (value_start <= value_end) {
                size_t part_end = attrs.find(',', value_start);
                if (part_end == std::string::npos || part_end > value_end) {
                    part_end = value_end;
                }
                if (part_end > value_start) {
                    values.push_back(attrs.substr(value_start, part_end - value_start));
                }
                if (part_end == value_end) {
                    break;
                }
                value_start = part_end + 1;
            }
            break;
        }

        pos = (value_end < attrs.size()) ? value_end + 1 : attrs.size();
    }
    return values;
}

static void print_qc_row(std::ostream& out,
                         const char* severity,
                         const char* code,
                         int line_idx,
                         const std::string& id,
                         const std::string& message) {
    out << severity << '\t'
        << code << '\t'
        << line_idx << '\t'
        << id << '\t'
        << message << '\n';
}

static void print_summary_tsv(std::ostream& out, const std::vector<SummaryRow>& rows) {
    out << "query_id\tmatched_id\tmatched_by\tseqid\tstart\tend\tstrand\ttype\tparent_id\t"
        << "child_count\ttranscript_count\texon_count\tcds_length\tstatus\n";
    for (const auto& row : rows) {
        out << row.query_id << '\t'
            << row.matched_id << '\t'
            << row.matched_by << '\t'
            << row.seqid << '\t'
            << row.start << '\t'
            << row.end << '\t'
            << row.strand << '\t'
            << row.type << '\t'
            << row.parent_id << '\t'
            << row.child_count << '\t'
            << row.transcript_count << '\t'
            << row.exon_count << '\t'
            << row.cds_length << '\t'
            << row.status << '\n';
    }
}

static void print_summary_json(std::ostream& out, const std::vector<SummaryRow>& rows) {
    out << "[\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        out << "  {"
            << "\"query_id\":\"" << json_escape(row.query_id) << "\","
            << "\"matched_id\":\"" << json_escape(row.matched_id) << "\","
            << "\"matched_by\":\"" << json_escape(row.matched_by) << "\","
            << "\"seqid\":\"" << json_escape(row.seqid) << "\","
            << "\"start\":" << row.start << ','
            << "\"end\":" << row.end << ','
            << "\"strand\":\"" << row.strand << "\","
            << "\"type\":\"" << json_escape(row.type) << "\","
            << "\"parent_id\":\"" << json_escape(row.parent_id) << "\","
            << "\"child_count\":" << row.child_count << ','
            << "\"transcript_count\":" << row.transcript_count << ','
            << "\"exon_count\":" << row.exon_count << ','
            << "\"cds_length\":" << row.cds_length << ','
            << "\"status\":\"" << row.status << "\""
            << "}";
        if (i + 1 < rows.size()) {
            out << ',';
        }
        out << '\n';
    }
    out << "]\n";
}

static int run_query(int argc, char* argv[], const char* prog) {
    if (argc < 2) {
        query_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::vector<std::string> ids;
    std::string name;
    std::string id_list_file;
    std::string region_str;
    std::string feature_type;
    std::string summary_format;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    bool include_children = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            ids.push_back(*value);
        } else if (arg == "--name") {
            auto value = require_value("--name");
            if (!value) return 1;
            name = *value;
        } else if (arg == "--id-list") {
            auto value = require_value("--id-list");
            if (!value) return 1;
            id_list_file = *value;
        } else if (arg == "--region") {
            auto value = require_value("--region");
            if (!value) return 1;
            region_str = *value;
        } else if (arg == "--type") {
            auto value = require_value("--type");
            if (!value) return 1;
            feature_type = *value;
        } else if (arg == "--attr") {
            auto value = require_value("--attr");
            if (!value) return 1;
            const auto equal_pos = value->find('=');
            if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value->size()) {
                std::cerr << "Error: --attr expects KEY=VALUE\n";
                return 1;
            }
            attr_filters.emplace_back(value->substr(0, equal_pos), value->substr(equal_pos + 1));
        } else if (arg == "--summary-format") {
            auto value = require_value("--summary-format");
            if (!value) return 1;
            summary_format = *value;
            if (summary_format != "tsv" && summary_format != "json") {
                std::cerr << "Error: --summary-format expects tsv or json\n";
                return 1;
            }
        } else if (arg == "--include-children") {
            include_children = true;
        } else if (arg == "-h" || arg == "--help") {
            query_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown query option " << arg << '\n';
            query_usage(prog);
            return 1;
        }
    }

    if (!id_list_file.empty()) {
        std::ifstream in{id_list_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) {
                ids.push_back(line);
            }
        }
    }

    gffsub::AnnotationIndex index = gffsub::AnnotationIndex::from_gff3(input_file);
    GffData result;
    std::unordered_set<int> seen;
    std::vector<SummaryRow> summary_rows;

    auto add_summary = [&](const std::string& query_id, const std::string& matched_by, const GffRecord& rec) {
        if (!summary_format.empty()) {
            summary_rows.push_back(make_summary_row(index, query_id, matched_by, rec));
        }
    };

    auto add_match = [&](const GffRecord& rec, const std::string& query_id, const std::string& matched_by) {
        if (feature_type.empty() || rec.type == feature_type) {
            if (append_unique(result, seen, rec)) {
                add_summary(query_id, matched_by, rec);
            }
        }
        if (include_children && rec.id) {
            for (const auto& child : index.descendants_of(*rec.id)) {
                if (feature_type.empty() || child.type == feature_type) {
                    if (append_unique(result, seen, child)) {
                        add_summary(query_id, "child", child);
                    }
                }
            }
        }
    };

    for (const auto& id : ids) {
        const auto rec = index.find_by_id(id);
        if (rec) {
            add_match(*rec, id, "ID");
        } else if (!summary_format.empty()) {
            summary_rows.push_back(make_not_found_row(id, "ID"));
        }
    }

    if (!name.empty()) {
        const auto rec = index.find_gene(name);
        if (rec) {
            add_match(*rec, name, infer_gene_match_key(index, name, *rec));
        } else if (!summary_format.empty()) {
            summary_rows.push_back(make_not_found_row(name, "name"));
        }
    }

    if (!region_str.empty()) {
        const auto region = parse_region(region_str);
        if (!region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        for (const auto& rec : index.overlap(region->seqid, region->start, region->end)) {
            add_match(rec, region_str, "region");
        }
    }

    for (const auto& [key, value] : attr_filters) {
        bool matched = false;
        for (const auto& rec : index.with_attribute(key, value)) {
            matched = true;
            add_match(rec, key + "=" + value, key);
        }
        if (!matched && !summary_format.empty()) {
            summary_rows.push_back(make_not_found_row(key + "=" + value, key));
        }
    }

    std::sort(result.records.begin(), result.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });
    if (summary_format == "tsv") {
        print_summary_tsv(std::cout, summary_rows);
    } else if (summary_format == "json") {
        print_summary_json(std::cout, summary_rows);
    } else {
        print_gff3(std::cout, result);
    }
    return 0;
}

static int run_window(int argc, char* argv[], const char* prog) {
    if (argc < 2) {
        window_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::string id;
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            id = *value;
        } else if (arg == "--upstream") {
            auto value = require_value("--upstream");
            if (!value) return 1;
            upstream = std::stoll(*value);
            if (upstream < 0) {
                std::cerr << "Error: --upstream must be non-negative\n";
                return 1;
            }
        } else if (arg == "--downstream") {
            auto value = require_value("--downstream");
            if (!value) return 1;
            downstream = std::stoll(*value);
            if (downstream < 0) {
                std::cerr << "Error: --downstream must be non-negative\n";
                return 1;
            }
        } else if (arg == "--strand-aware") {
            strand_aware = true;
        } else if (arg == "-h" || arg == "--help") {
            window_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown window option " << arg << '\n';
            window_usage(prog);
            return 1;
        }
    }

    if (id.empty()) {
        std::cerr << "Error: window requires --id\n";
        return 1;
    }

    const auto index = gffsub::AnnotationIndex::from_gff3(input_file);
    auto target = index.find_by_id(id);
    if (!target) {
        target = index.find_gene(id);
    }
    if (!target) {
        std::cerr << "Error: cannot find " << id << '\n';
        return 1;
    }

    const auto region = window_region(*target, upstream, downstream, strand_aware);
    GffData result;
    for (const auto& rec : index.overlap(region.seqid, region.start, region.end)) {
        result.append(rec);
    }
    print_gff3(std::cout, result);
    return 0;
}

static int run_qc(int argc, char* argv[], const char* prog) {
    if (argc != 2) {
        qc_usage(prog);
        return 1;
    }

    GffData data;
    IdIndex idx;
    const std::string input_file = argv[1];
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    std::unordered_map<std::string, const GffRecord*> by_id;
    std::unordered_map<std::string, int> id_counts;
    for (const auto& rec : data.records) {
        if (rec.id) {
            ++id_counts[*rec.id];
            by_id.emplace(*rec.id, &rec);
        }
    }

    std::cout << "severity\tcode\tline_idx\tid\tmessage\n";

    for (const auto& [id, count] : id_counts) {
        if (count > 1) {
            print_qc_row(std::cout, "error", "duplicate_id", -1, id, "ID appears more than once");
        }
    }

    for (const auto& rec : data.records) {
        const std::string id = record_id(rec);
        if (rec.start > rec.end) {
            print_qc_row(std::cout, "error", "invalid_range", rec.line_idx, id, "start is greater than end");
        }

        for (const auto& parent_id : attribute_values(rec.attr_raw, "Parent")) {
            const auto parent_it = by_id.find(parent_id);
            if (parent_it == by_id.end()) {
                print_qc_row(std::cout, "error", "missing_parent", rec.line_idx, id,
                             "Parent " + parent_id + " was not found");
                continue;
            }

            const auto& parent = *parent_it->second;
            if (rec.seqid != parent.seqid || rec.start < parent.start || rec.end > parent.end) {
                print_qc_row(std::cout, "warning", "child_outside_parent", rec.line_idx, id,
                             "child is outside Parent " + parent_id);
            }
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "query") {
        return run_query(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "window") {
        return run_window(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "qc") {
        return run_qc(argc - 1, argv + 1, argv[0]);
    }

    std::string region_str;
    std::string bed_file;
    std::string feature;
    bool do_longest = false;
    size_t num_threads = 6;
    std::string output_format = "gff3";
    std::string output_file;

    static struct option long_options[] = {
        {"region",        required_argument, nullptr, 'r'},
        {"bed",           required_argument, nullptr, 'b'},
        {"feature",       required_argument, nullptr, 'f'},
        {"longest",       no_argument,       nullptr, 'L'},
        {"threads",       required_argument, nullptr, '@'},
        {"output-format", required_argument, nullptr, 't'},
        {"output",        required_argument, nullptr, 'o'},
        {"help",          no_argument,       nullptr, 'h'},
        {nullptr,        0,                 nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "r:b:f:L@:t:o:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'r': region_str = optarg; break;
            case 'b': bed_file = optarg; break;
            case 'f': feature = optarg; break;
            case 'L': do_longest = true; break;
            case '@': {
                size_t t = std::stoul(optarg);
                if (t == 0) t = 1;
                if (t > 256) t = 256; // cap to prevent over-subscription
                num_threads = t;
                break;
            }
            case 't': output_format = optarg; break;
            case 'o': output_file = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }

    // Validate output format
    OutputFormat fmt = OutputFormat::GFF3;
    if (output_format == "gtf") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf2") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf3") {
        fmt = OutputFormat::GTF3;
    } else if (output_format == "bed") {
        fmt = OutputFormat::BED;
    } else if (output_format == "gff3") {
        fmt = OutputFormat::GFF3;
    } else {
        std::cerr << "Error: unknown output format " << output_format << '\n';
        std::cerr << "Supported formats: gff3, gtf2, gtf3, bed\n";
        return 1;
    }

    std::string input_file = argv[optind];

    GffData data;
    IdIndex idx;

    // Parse input file
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    Region region{"", 0, 0};
    std::optional<Region> parsed_region;

    // Apply region filters
    if (!region_str.empty()) {
        parsed_region = parse_region(region_str);
        if (!parsed_region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        region = *parsed_region;
        filter_by_region(data, region);
    }

    if (!bed_file.empty()) {
        filter_by_regions_from_file(data, bed_file);
    }

    // Apply feature filters
    if (do_longest) {
        filter_longest(data, idx, feature, num_threads);
    } else if (!feature.empty()) {
        filter_by_feature(data, feature);
    }

    // Output
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!output_file.empty()) {
        out_file.open(output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: cannot open " << output_file << '\n';
            return 1;
        }
        out = &out_file;
    }

    switch (fmt) {
        case OutputFormat::GFF3: print_gff3(*out, data); break;
        case OutputFormat::GTF2:
        case OutputFormat::GTF3: print_gtf(*out, data, fmt); break;
        case OutputFormat::BED:  print_bed(*out, data); break;
    }

    return 0;
}
