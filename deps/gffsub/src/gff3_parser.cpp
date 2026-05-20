#include "gff3.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>

namespace gffsub {

static std::optional<std::string> extract_attr_value(std::string_view attrs, std::string_view key) {
    size_t pos = 0;
    while (pos < attrs.size()) {
        size_t key_end = attrs.find('=', pos);
        if (key_end == std::string_view::npos) {
            size_t next = attrs.find(';', pos);
            pos = (next == std::string_view::npos) ? attrs.size() : next + 1;
            continue;
        }

        std::string_view found_key = attrs.substr(pos, key_end - pos);
        size_t val_start = key_end + 1;
        size_t val_end = attrs.find(';', val_start);
        if (val_end == std::string_view::npos) val_end = attrs.size();

        if (found_key == key) {
            return std::string(attrs.substr(val_start, val_end - val_start));
        }
        pos = (val_end < attrs.size()) ? val_end + 1 : attrs.size();
    }
    return std::nullopt;
}

static std::optional<std::string> extract_quoted_value(const std::string& attrs, const char* key) {
    size_t pos = attrs.find(key);
    if (pos == std::string::npos) return std::nullopt;
    size_t q1 = attrs.find('"', pos);
    if (q1 == std::string::npos) return std::nullopt;
    size_t q2 = attrs.find('"', q1 + 1);
    if (q2 == std::string::npos) return std::nullopt;
    return attrs.substr(q1 + 1, q2 - q1 - 1);
}

static std::vector<std::string> split_line(const std::string& line, char delimiter) {
    std::vector<std::string> cols;
    size_t start = 0;
    while (true) {
        auto pos = line.find(delimiter, start);
        if (pos == std::string::npos) {
            cols.emplace_back(line.substr(start));
            break;
        }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

int parse_file(const std::string& filename, GffData& data, IdIndex& idx, InputFormat format) {
    std::ifstream file(filename);
    if (!file.is_open()) return -1;

    std::string line;
    int line_num = 0;
    bool in_fasta = false;

    while (std::getline(file, line)) {
        line_num++;
        if (in_fasta) continue;

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty() || line[0] == '#') continue;

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) continue;

            rec.seqid = cols[0];
            rec.source = cols[1];
            rec.type = cols[2];
            rec.start = std::stoll(cols[3]);
            rec.end = std::stoll(cols[4]);
            rec.score = (cols[5] == ".") ? std::nullopt : std::optional(std::stod(cols[5]));
            rec.strand = cols[6].empty() ? '.' : cols[6][0];
            rec.phase = cols[7].empty() ? '.' : cols[7][0];
            rec.attr_raw = cols[8];

            rec.id = extract_attr_value(cols[8], "ID");
            rec.parent_id = extract_attr_value(cols[8], "Parent");
            rec.gene_id = extract_attr_value(cols[8], "gene_id");
            rec.transcript_id = extract_attr_value(cols[8], "transcript_id");

            if (format == InputFormat::GTF) {
                if (!rec.gene_id) {
                    rec.gene_id = extract_quoted_value(cols[8], "gene_id");
                }
                if (!rec.transcript_id) {
                    rec.transcript_id = extract_quoted_value(cols[8], "transcript_id");
                }
            }
        } else if (format == InputFormat::BED) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            rec.seqid = cols[0];
            rec.start = std::stoll(cols[1]) + 1;
            rec.end = std::stoll(cols[2]);
            rec.source = "gffsub";
            rec.type = "region";
            rec.score = cols.size() > 4 ? std::optional(std::stod(cols[4])) : std::nullopt;
            rec.strand = cols.size() > 5 ? (cols[5][0]) : '.';
            rec.phase = '.';
            rec.id = cols.size() > 3 ? std::optional(cols[3]) : std::nullopt;
            rec.parent_id = std::nullopt;
            rec.gene_id = std::nullopt;
            rec.transcript_id = std::nullopt;
            rec.attr_raw = rec.id ? "ID=" + *rec.id : "";
        }

        if (rec.id) {
            idx.add(*rec.id, rec.line_idx);
        }

        data.append(rec);
    }

    return 0;
}

std::optional<Region> parse_region(std::string_view region_str) {
    size_t colon = region_str.find(':');
    if (colon == std::string_view::npos) return std::nullopt;

    std::string seqid(region_str.substr(0, colon));
    auto range_part = region_str.substr(colon + 1);

    size_t dash = range_part.find('-');
    if (dash == std::string_view::npos) return std::nullopt;

    int64_t start = std::stoll(std::string(range_part.substr(0, dash)));
    int64_t end = std::stoll(std::string(range_part.substr(dash + 1)));

    return Region{seqid, start, end};
}

BedRegion to_bed_region(const GffRecord& rec) {
    return BedRegion{rec.seqid, rec.start - 1, rec.end};
}

Region from_bed_region(const BedRegion& region) {
    return Region{region.seqid, region.start + 1, region.end};
}

Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware) {
    int64_t left_extension = upstream;
    int64_t right_extension = downstream;
    if (strand_aware && rec.strand == '-') {
        left_extension = downstream;
        right_extension = upstream;
    }

    int64_t start = rec.start - left_extension;
    if (start < 1) {
        start = 1;
    }
    return Region{rec.seqid, start, rec.end + right_extension};
}

void filter_by_region(GffData& data, const Region& region) {
    for (auto& rec : data) {
        if (rec.seqid != region.seqid || rec.end < region.start || rec.start > region.end) {
            rec.kept = false;
        }
    }
}

static std::vector<Region> load_regions(const std::string& filename, bool is_bed) {
    std::vector<Region> regions;
    if (is_bed) {
        std::ifstream file(filename);
        if (!file.is_open()) return regions;

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            Region r;
            r.seqid = cols[0];
            r.start = std::stoll(cols[1]) + 1;
            r.end = std::stoll(cols[2]);
            regions.push_back(r);
        }
    }
    return regions;
}

void filter_by_regions_from_file(GffData& data, const std::string& bed_file) {
    auto regions = load_regions(bed_file, true);
    for (auto& rec : data) {
        bool in_region = false;
        for (const auto& r : regions) {
            if (rec.seqid == r.seqid && rec.end >= r.start && rec.start <= r.end) {
                in_region = true;
                break;
            }
        }
        if (!in_region) rec.kept = false;
    }
}

void filter_by_feature(GffData& data, std::string_view feature_type) {
    for (auto& rec : data) {
        if (rec.kept && rec.type != feature_type) {
            rec.kept = false;
        }
    }
}

}  // namespace gffsub
