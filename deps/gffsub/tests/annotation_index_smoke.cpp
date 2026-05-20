#include "gff3.hpp"

#include <fstream>
#include <iostream>
#include <string>

static bool write_test_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=GeneOne;gene_id=G1;locus_tag=Locus1;Alias=Alpha,Beta;Dbxref=GeneID:12345\n"
        << "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
        << "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx2;Parent=gene1\n"
        << "chr1\t.\texon\t10\t20\t.\t+\t.\tID=exon_shared;Parent=tx1,tx2\n"
        << "chr1\t.\tCDS\t30\t50\t.\t+\t0\tID=cds1;Parent=tx1\n"
        << "chr1\t.\tgene\t200\t300\t.\t-\t.\tID=gene2;Name=GeneTwo\n";
    return true;
}

static int check_self_contained_annotation() {
    const std::string path{"annotation_index_smoke.gff3"};
    if (!write_test_annotation(path)) {
        std::cerr << "cannot write self-contained test annotation\n";
        return 1;
    }

    const auto index = gffsub::AnnotationIndex::from_gff3(path);

    const auto by_id = index.find_by_id("gene1");
    if (!by_id || by_id->type != "gene" || by_id->seqid != "chr1") {
        std::cerr << "find_by_id failed for gene1\n";
        return 1;
    }

    const auto bed_region = gffsub::to_bed_region(*by_id);
    if (bed_region.seqid != "chr1" || bed_region.start != 0 || bed_region.end != 100) {
        std::cerr << "to_bed_region failed for gene1\n";
        return 1;
    }

    const auto gff_region = gffsub::from_bed_region(bed_region);
    if (gff_region.seqid != "chr1" || gff_region.start != 1 || gff_region.end != 100) {
        std::cerr << "from_bed_region failed for gene1\n";
        return 1;
    }

    const auto plus_window = gffsub::window_region(*by_id, 10, 20, true);
    if (plus_window.seqid != "chr1" || plus_window.start != 1 || plus_window.end != 120) {
        std::cerr << "window_region failed for plus-strand gene1\n";
        return 1;
    }

    for (const std::string query : {"GeneOne", "G1", "Locus1", "Alpha", "Beta", "GeneID:12345"}) {
        const auto gene = index.find_gene(query);
        if (!gene || gene->id != "gene1") {
            std::cerr << "find_gene failed for " << query << '\n';
            return 1;
        }
    }

    const auto children = index.children_of("gene1");
    if (children.size() != 2) {
        std::cerr << "children_of failed for gene1\n";
        return 1;
    }

    const auto parents = index.parents_of("exon_shared");
    if (parents.size() != 2) {
        std::cerr << "multi-parent parents_of failed\n";
        return 1;
    }

    const auto descendants = index.descendants_of("gene1");
    if (descendants.size() != 5) {
        std::cerr << "descendants_of failed for gene1\n";
        return 1;
    }

    const auto overlapping = index.overlap("chr1", 15, 35);
    if (overlapping.size() != 5) {
        std::cerr << "overlap failed for chr1:15-35\n";
        return 1;
    }

    const auto attr_matches = index.with_attribute("Alias", "Beta");
    if (attr_matches.size() != 1 || attr_matches.front().id != "gene1") {
        std::cerr << "with_attribute failed for Alias=Beta\n";
        return 1;
    }

    const auto nearest = index.nearest_gene("chr1", 120, 150);
    if (!nearest || nearest->id != "gene1") {
        std::cerr << "nearest_gene failed for chr1:120-150\n";
        return 1;
    }

    const auto model = index.gene_model("exon_shared");
    if (!model || model->gene.id != "gene1" || model->records.size() != 5) {
        std::cerr << "gene_model failed for exon_shared\n";
        return 1;
    }
    if (model->records.front().id != "gene1" || model->records.back().id != "cds1") {
        std::cerr << "gene_model record order failed\n";
        return 1;
    }

    const auto minus_gene = index.find_by_id("gene2");
    if (!minus_gene) {
        std::cerr << "find_by_id failed for gene2\n";
        return 1;
    }
    const auto minus_window = gffsub::window_region(*minus_gene, 10, 20, true);
    if (minus_window.seqid != "chr1" || minus_window.start != 180 || minus_window.end != 310) {
        std::cerr << "window_region failed for minus-strand gene2\n";
        return 1;
    }

    return 0;
}

static int check_soybean_annotation(const std::string& path) {
    const std::string gene_id{"SoyL04_01G000000"};
    const auto index = gffsub::AnnotationIndex::from_gff3(path);

    const auto gene = index.find_by_id(gene_id);
    if (!gene || gene->type != "gene" || gene->seqid != "Chr01") {
        std::cerr << "find_by_id failed for " << gene_id << '\n';
        return 1;
    }

    const auto gene_by_lookup = index.find_gene(gene_id);
    if (!gene_by_lookup || gene_by_lookup->id != gene->id) {
        std::cerr << "find_gene failed for " << gene_id << '\n';
        return 1;
    }

    const auto children = index.children_of(gene_id);
    if (children.size() != 1 || children.front().type != "mRNA") {
        std::cerr << "children_of failed for " << gene_id << '\n';
        return 1;
    }

    const auto parents = index.parents_of("SoyL04_01G000000.m1");
    if (parents.size() != 1 || parents.front().id != gene->id) {
        std::cerr << "parents_of failed for SoyL04_01G000000.m1\n";
        return 1;
    }

    const auto descendants = index.descendants_of(gene_id);
    if (descendants.size() != 7) {
        std::cerr << "descendants_of failed for " << gene_id << '\n';
        return 1;
    }

    const auto overlapping = index.overlap("Chr01", 1, 35000);
    if (overlapping.size() != 8 || overlapping.front().id != gene->id) {
        std::cerr << "overlap failed for Chr01:1-35000\n";
        return 1;
    }

    const auto attr_matches = index.with_attribute("ID", gene_id);
    if (attr_matches.size() != 1 || attr_matches.front().id != gene->id) {
        std::cerr << "with_attribute failed for ID=" << gene_id << '\n';
        return 1;
    }

    const auto nearest = index.nearest_gene("Chr01", 35001, 36000);
    if (!nearest || nearest->id != gene->id) {
        std::cerr << "nearest_gene failed for Chr01:35001-36000\n";
        return 1;
    }

    const auto model = index.gene_model("SoyL04_01G000000.m1");
    if (!model || model->gene.id != gene->id || model->records.size() != 8) {
        std::cerr << "gene_model failed for SoyL04_01G000000.m1\n";
        return 1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 2) {
        std::cerr << "usage: annotation_index_smoke [annotation.gff3]\n";
        return 2;
    }

    if (argc == 2) {
        if (check_soybean_annotation(argv[1]) != 0) {
            return 1;
        }
    } else {
        const int self_contained_status = check_self_contained_annotation();
        if (self_contained_status != 0) {
            return self_contained_status;
        }
    }

    std::cout << "annotation_index_smoke OK\n";
    return 0;
}
