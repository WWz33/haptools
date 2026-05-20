#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "gff3.hpp"

namespace haplokit {

struct GeneAnnotation {
    std::string mode;  // "overlap", "nearest", "gene", "none"
    std::string id;
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::optional<int64_t> distance;
};

class GffAnnotator {
public:
    bool load(const std::string& gff_path);
    GeneAnnotation annotate(const std::string& chrom, int64_t start, int64_t end) const;
    std::optional<GeneAnnotation> find_gene(const std::string& gene_id) const;

private:
    std::optional<gffsub::AnnotationIndex> index_;
};

}  // namespace haplokit
