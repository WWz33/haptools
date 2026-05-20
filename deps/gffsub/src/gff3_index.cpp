#include "gff3.hpp"

#include <deque>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace gffsub {

static std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs) {
    std::unordered_map<std::string, std::vector<std::string>> parsed;
    size_t pos = 0;
    while (pos < attrs.size()) {
        const size_t key_end = attrs.find('=', pos);
        if (key_end == std::string_view::npos) {
            break;
        }

        const std::string key{attrs.substr(pos, key_end - pos)};
        const size_t value_start = key_end + 1;
        size_t value_end = attrs.find(';', value_start);
        if (value_end == std::string_view::npos) {
            value_end = attrs.size();
        }

        const std::string value{attrs.substr(value_start, value_end - value_start)};
        if (!key.empty() && !value.empty()) {
            size_t part_start = 0;
            while (part_start <= value.size()) {
                size_t part_end = value.find(',', part_start);
                if (part_end == std::string::npos) {
                    part_end = value.size();
                }
                const std::string part = value.substr(part_start, part_end - part_start);
                if (!part.empty()) {
                    parsed[key].push_back(part);
                }
                if (part_end == value.size()) {
                    break;
                }
                part_start = part_end + 1;
            }
        }

        pos = (value_end < attrs.size()) ? value_end + 1 : attrs.size();
    }
    return parsed;
}

AnnotationIndex AnnotationIndex::from_gff3(const std::string& path) {
    GffData data;
    IdIndex idx;
    if (parse_file(path, data, idx, InputFormat::GFF3) != 0) {
        throw std::runtime_error("cannot parse GFF3 file: " + path);
    }
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex::AnnotationIndex(GffData data) : data_(std::move(data)) {
    for (int i = 0; i < static_cast<int>(data_.records.size()); ++i) {
        const auto& rec = data_.records[i];
        if (rec.id) {
            id_to_record_.emplace(*rec.id, i);
        }

        if (rec.type != "gene") {
            continue;
        }

        auto add_gene_key = [&](const std::string& key) {
            if (!key.empty()) {
                gene_lookup_[key].push_back(i);
            }
        };

        if (rec.id) {
            add_gene_key(*rec.id);
        }
        if (rec.gene_id) {
            add_gene_key(*rec.gene_id);
        }

        const auto attrs = parse_attributes(rec.attr_raw);
        for (const char* key : {"Name", "gene_id", "locus_tag", "Alias", "Dbxref"}) {
            const auto it = attrs.find(key);
            if (it == attrs.end()) {
                continue;
            }
            for (const auto& value : it->second) {
                add_gene_key(value);
            }
        }
    }

    for (int i = 0; i < static_cast<int>(data_.records.size()); ++i) {
        const auto& rec = data_.records[i];
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto parent_it = attrs.find("Parent");
        if (parent_it == attrs.end()) {
            continue;
        }

        for (const auto& parent_id : parent_it->second) {
            children_by_parent_id_[parent_id].push_back(i);
            const auto parent_record = id_to_record_.find(parent_id);
            if (rec.id && parent_record != id_to_record_.end()) {
                parents_by_child_id_[*rec.id].push_back(parent_record->second);
            }
        }
    }
}

std::optional<GffRecord> AnnotationIndex::find_by_id(std::string_view id) const {
    const auto it = id_to_record_.find(std::string{id});
    if (it == id_to_record_.end()) {
        return std::nullopt;
    }
    return data_.records[it->second];
}

std::optional<GffRecord> AnnotationIndex::find_gene(std::string_view id) const {
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end() || it->second.empty()) {
        return std::nullopt;
    }
    return data_.records[it->second.front()];
}

std::vector<GffRecord> AnnotationIndex::parents_of(std::string_view id) const {
    std::vector<GffRecord> parents;
    const auto it = parents_by_child_id_.find(std::string{id});
    if (it == parents_by_child_id_.end()) {
        return parents;
    }
    parents.reserve(it->second.size());
    for (const int idx : it->second) {
        parents.push_back(data_.records[idx]);
    }
    return parents;
}

std::vector<GffRecord> AnnotationIndex::children_of(std::string_view parent_id) const {
    std::vector<GffRecord> children;
    const auto it = children_by_parent_id_.find(std::string{parent_id});
    if (it == children_by_parent_id_.end()) {
        return children;
    }
    children.reserve(it->second.size());
    for (const int idx : it->second) {
        children.push_back(data_.records[idx]);
    }
    return children;
}

std::vector<GffRecord> AnnotationIndex::descendants_of(std::string_view parent_id) const {
    std::vector<GffRecord> descendants;
    std::deque<std::string> pending{std::string{parent_id}};
    std::unordered_set<std::string> visited;

    while (!pending.empty()) {
        const auto current = pending.front();
        pending.pop_front();
        if (!visited.insert(current).second) {
            continue;
        }

        const auto child_it = children_by_parent_id_.find(current);
        if (child_it == children_by_parent_id_.end()) {
            continue;
        }

        for (const int child_idx : child_it->second) {
            const auto& child = data_.records[child_idx];
            descendants.push_back(child);
            if (child.id) {
                pending.push_back(*child.id);
            }
        }
    }

    return descendants;
}

std::vector<GffRecord> AnnotationIndex::overlap(std::string_view seqid, int64_t start, int64_t end) const {
    std::vector<GffRecord> matches;
    for (const auto& rec : data_.records) {
        if (rec.seqid == seqid && rec.end >= start && rec.start <= end) {
            matches.push_back(rec);
        }
    }
    return matches;
}

std::optional<GffRecord> AnnotationIndex::nearest_gene(std::string_view seqid, int64_t start, int64_t end) const {
    std::optional<GffRecord> nearest;
    int64_t nearest_distance = std::numeric_limits<int64_t>::max();

    for (const auto& rec : data_.records) {
        if (rec.type != "gene" || rec.seqid != seqid) {
            continue;
        }

        int64_t distance = 0;
        if (end < rec.start) {
            distance = rec.start - end;
        } else if (start > rec.end) {
            distance = start - rec.end;
        }

        if (distance < nearest_distance) {
            nearest_distance = distance;
            nearest = rec;
        }
    }

    return nearest;
}

std::vector<GffRecord> AnnotationIndex::with_attribute(std::string_view key, std::string_view value) const {
    std::vector<GffRecord> matches;
    for (const auto& rec : data_.records) {
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto it = attrs.find(std::string{key});
        if (it == attrs.end()) {
            continue;
        }
        for (const auto& attr_value : it->second) {
            if (attr_value == value) {
                matches.push_back(rec);
                break;
            }
        }
    }
    return matches;
}

std::optional<GeneModel> AnnotationIndex::gene_model(std::string_view id) const {
    std::optional<GffRecord> gene = find_gene(id);
    if (!gene) {
        auto rec = find_by_id(id);
        while (rec && rec->type != "gene" && rec->id) {
            const auto parents = parents_of(*rec->id);
            if (parents.empty()) {
                rec = std::nullopt;
            } else {
                rec = parents.front();
            }
        }
        if (rec && rec->type == "gene") {
            gene = rec;
        }
    }

    if (!gene || !gene->id) {
        return std::nullopt;
    }

    GeneModel model{*gene, {}};
    model.records.push_back(*gene);

    std::unordered_set<int> seen{gene->line_idx};
    for (const auto& rec : descendants_of(*gene->id)) {
        if (seen.insert(rec.line_idx).second) {
            model.records.push_back(rec);
        }
    }

    std::sort(model.records.begin(), model.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });

    return model;
}

}  // namespace gffsub
