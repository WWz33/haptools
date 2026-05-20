#include "gff_annotator.h"

#include <exception>

namespace haplokit {

namespace {

GeneAnnotation to_gene_annotation(const gffsub::GffRecord& rec, const std::string& mode,
                                  std::optional<int64_t> distance = std::nullopt) {
    GeneAnnotation ann;
    ann.mode = mode;
    ann.id = rec.id.value_or(rec.gene_id.value_or(""));
    ann.seqid = rec.seqid;
    ann.start = rec.start;
    ann.end = rec.end;
    ann.strand = rec.strand;
    ann.distance = distance;
    return ann;
}

int64_t distance_to_record(const gffsub::GffRecord& rec, int64_t start, int64_t end) {
    if (end < rec.start) {
        return rec.start - end;
    }
    if (start > rec.end) {
        return start - rec.end;
    }
    return 0;
}

}  // namespace

bool GffAnnotator::load(const std::string& gff_path) {
    try {
        index_ = gffsub::AnnotationIndex::from_gff3(gff_path);
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

GeneAnnotation GffAnnotator::annotate(const std::string& chrom, int64_t start, int64_t end) const {
    GeneAnnotation ann;
    if (!index_) {
        ann.mode = "none";
        return ann;
    }

    for (const auto& rec : index_->overlap(chrom, start, end)) {
        if (rec.type == "gene") {
            return to_gene_annotation(rec, "overlap");
        }
    }

    const auto nearest = index_->nearest_gene(chrom, start, end);
    if (nearest) {
        return to_gene_annotation(*nearest, "nearest", distance_to_record(*nearest, start, end));
    }

    ann.mode = "none";
    return ann;
}

std::optional<GeneAnnotation> GffAnnotator::find_gene(const std::string& gene_id) const {
    if (!index_) {
        return std::nullopt;
    }
    const auto gene = index_->find_gene(gene_id);
    if (!gene) {
        return std::nullopt;
    }
    return to_gene_annotation(*gene, "gene");
}

std::optional<GeneAnnotation> GffAnnotator::find_gene_window(
    const std::string& gene_id, int64_t upstream, int64_t downstream, bool strand_aware) const {
    if (!index_) {
        return std::nullopt;
    }
    const auto gene = index_->find_gene(gene_id);
    if (!gene) {
        return std::nullopt;
    }

    auto ann = to_gene_annotation(*gene, "gene");
    const auto region = gffsub::window_region(*gene, upstream, downstream, strand_aware);
    ann.seqid = region.seqid;
    ann.start = region.start;
    ann.end = region.end;
    return ann;
}

}  // namespace haplokit
