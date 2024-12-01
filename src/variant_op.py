import numpy as np
import pysam


class SimpleVariant:
    # # Variant types: DEL, INS, INV, tDUP, itDUP, dDUP, idDUP, SNP
    def __init__(self, chrom, start, end, type, alt_bases, supp_reads, supp_cigar):

        self.ref_chrom = chrom
        self.ref_start = start
        self.ref_end = end

        self.type = type
        self.alt_bases = alt_bases

        self.supp_reads = supp_reads
        self.supp_cigar = supp_cigar

        self.hybrid_length = supp_cigar.hybrid_length
        self.read_length = self.calculate_read_length()
        self.ref_length = self.ref_end - self.ref_start

        self.dup_source_chrom = None
        self.dup_source_start = None
        self.dup_source_end = None
        self.dup_source_strand = None

        self.bnd_source_chrom = None
        self.bnd_source_start = None
        self.bnd_source_end = None
        self.bnd_source_strand = None

        self.pseudo_vaf = 0
        self.genuine_vaf = 0

    def set_pseudo_vaf(self, vaf):
        self.pseudo_vaf = vaf
        self.genuine_vaf = self.pseudo_vaf

    def calculate_read_length(self):
        read_length = 0

        if self.supp_cigar.op != "D":
            read_length += self.supp_cigar.hybrid_length

        return read_length

    def set_duplicated_source(self, chrom, start, end, strand):
        self.dup_source_chrom = chrom
        self.dup_source_start = start
        self.dup_source_end = end
        self.dup_source_strand = strand

    def set_bnd_source(self, chrom, start, end, strand):
        self.bnd_source_chrom = chrom
        self.bnd_source_start = start
        self.bnd_source_end = end
        self.bnd_source_strand = strand

    def to_string(self):
        return "{}-{}-{},length={},type={},supp={},seq={},DUP={}-{}-{},BND={}-{}-{}".format(self.ref_chrom, self.ref_start, self.ref_end, self.hybrid_length, self.type, len(self.supp_reads), "N", self.dup_source_chrom, self.dup_source_start, self.dup_source_end, self.bnd_source_chrom, self.bnd_source_start, self.bnd_source_end)


def detect_variant_from_detail_cigar(detail_cigar, supp_reads, options):
    """
    detect simple variant from detail cigar
    """

    # # each cigar represent a simple variant
    variant_chrom, variant_start, variant_end, variant_strand, variant_supp_reads = detail_cigar.ref_chrom, detail_cigar.ref_start, detail_cigar.ref_end, detail_cigar.strand, supp_reads

    if detail_cigar.op == "B":
        variant = SimpleVariant(variant_chrom, variant_start, variant_end, "BND", "N", variant_supp_reads, detail_cigar)
        variant.set_bnd_source(detail_cigar.bnd_ref_chrom, detail_cigar.bnd_ref_start, detail_cigar.bnd_ref_end, detail_cigar.bnd_ref_strand)

    elif detail_cigar.op == "D":
        variant = SimpleVariant(variant_chrom, variant_start, variant_end, "DEL", "N", variant_supp_reads, detail_cigar)

    elif detail_cigar.op == "V":
        variant = SimpleVariant(variant_chrom, variant_start, variant_end, "INV", "N", variant_supp_reads, detail_cigar)

    elif detail_cigar.op in ["I", "MappedI", "UnmappedI"]:
        if detail_cigar.op == "UnmappedI":
            variant = SimpleVariant(variant_chrom, variant_start, variant_end, "INS", "N", variant_supp_reads, detail_cigar)
        else:
            # # STEP: add DUP, including both tDUP and dDUP
            mapped_ref_chrom, mapped_ref_start, mapped_ref_end, mapped_strand = detail_cigar.secondary_ref_chrom, detail_cigar.secondary_ref_start, detail_cigar.secondary_ref_end, detail_cigar.secondary_ref_strand

            if mapped_ref_chrom != variant_chrom:
                # variant_type = "dDUP/BND" if mapped_strand == "+" else "idDUP/BND"
                variant_type = "dDUP" if mapped_strand == "+" else "idDUP"

            else:
                # # determine the variant type
                if abs(mapped_ref_start - variant_start) <= options.min_sv_size or abs(mapped_ref_end - variant_start) <= options.min_sv_size:
                    variant_type = "tDUP" if mapped_strand == "+" else "itDUP"
                else:
                    if abs(variant_start - mapped_ref_start) > options.max_sv_size:
                        # variant_type = "dDUP/BND" if mapped_strand == "+" else "idDUP/BND"
                        variant_type = "dDUP" if mapped_strand == "+" else "idDUP"

                    else:
                        variant_type = "dDUP" if mapped_strand == "+" else "idDUP"

            variant = SimpleVariant(variant_chrom, variant_start, variant_end, variant_type, "N", variant_supp_reads, detail_cigar)  # -1: convert 1-based to 0-based, +1: convert [) to []
            if "BND" in variant_type:
                variant.set_bnd_source(mapped_ref_chrom, mapped_ref_start, mapped_ref_end, mapped_strand)
            variant.set_duplicated_source(mapped_ref_chrom, mapped_ref_start, mapped_ref_end, mapped_strand)

    # elif detail_cigar.op in ["A", "T", "C", "G", "X"]:
    #
    #     variant = SimpleVariant(variant_chrom, variant_start, variant_end, "SNP", "N", variant_supp_reads, detail_cigar)

    else:
        variant = None
        print("[ERROR]: no this op", detail_cigar.op)
        exit()

    variant.set_pseudo_vaf(detail_cigar.pseudo_vaf)

    return variant


def detect_pseudo_vaf_for_detail_cigar_via_coverage_list(detail_cigar, supp_num, coverage_list_ref_start, coverage_list_ref_end, coverage_list, options):
    """
    detect pseudo vaf for detail cigar by calculating the read depth from bam
    """

    if not coverage_list_ref_start < detail_cigar.ref_start < coverage_list_ref_end:
        detect_pseudo_vaf_for_detail_cigar_via_fetch_bam(detail_cigar, supp_num, options.target_path, options.min_mapq)

        # unique_read_num = supp_num
        # # # STEP: calculate pseudo vaf
        # try:
        #     pseudo_vaf = round(supp_num / unique_read_num, 2)
        # except ZeroDivisionError:
        #     pseudo_vaf = 1.0
        #
        # # print(supp_num, unique_read_num, pseudo_vaf)
        # # print("not covered", detail_cigar.ref_start, supp_num, unique_read_num, pseudo_vaf)
        #
        # # # STEP: set
        # detail_cigar.set_coverage(unique_read_num)
        # detail_cigar.set_pseudo_vaf(min(1.0, pseudo_vaf))
        
    else:

        # unique_read_num = int(coverage_list[detail_cigar.ref_start - coverage_list_ref_start])
        coverage_list_len = len(coverage_list)
        unique_read_num = max([int(coverage_list[min(detail_cigar.ref_start + i - coverage_list_ref_start, coverage_list_len - 1)]) for i in range(-5, 5)])

        # # STEP: calculate pseudo vaf
        try:
            pseudo_vaf = round(supp_num / unique_read_num, 2)
        except ZeroDivisionError:
            pseudo_vaf = 1.0
        # print("covered", detail_cigar.ref_start, supp_num, unique_read_num, pseudo_vaf)
        # print(supp_num, unique_read_num, pseudo_vaf)

        # for i in range(-5, 5):
        #     print(int(coverage_list[detail_cigar.ref_start + i - coverage_list_ref_start]))

        # # STEP: set
        detail_cigar.set_coverage(unique_read_num)
        detail_cigar.set_pseudo_vaf(min(1.0, pseudo_vaf))


def detect_pseudo_vaf_for_detail_cigar_via_fetch_bam(detail_cigar, supp_num, bam_path, min_mapq):
    """
    detect pseudo vaf for detail cigar by calculating the read depth from bam
    """

    bam_file = pysam.AlignmentFile(bam_path)

    # # STEP: collect unique read number at both left and right bkp
    unique_read_list = []

    for bkp_pos in [detail_cigar.ref_start, detail_cigar.ref_end]:
        left_extend_bkp_pos = bkp_pos - 10
        right_extend_bkp_pos = bkp_pos + 10

        for align in bam_file.fetch(detail_cigar.ref_chrom, bkp_pos - 2, bkp_pos + 2):
            if align.qname not in unique_read_list and align.mapq >= min_mapq and align.reference_start < left_extend_bkp_pos and align.reference_end > right_extend_bkp_pos:
                unique_read_list.append(align.qname)

    unique_read_num = len(unique_read_list)

    # # STEP: calculate pseudo vaf
    try:
        pseudo_vaf = round(supp_num / unique_read_num, 2)
    except ZeroDivisionError:
        pseudo_vaf = 1.0

    # # STEP: set
    detail_cigar.set_coverage(unique_read_num)
    detail_cigar.set_pseudo_vaf(min(1.0, pseudo_vaf))


def detect_target_pseudo_variant_for_partition(partition, interval_start, interval_end, target_interval_coverage, options):
    """
    segment and detect SV from single read's non_match_cigar
    one target variant may have more than one base variants
    :return:
    """

    # # STEP: detect variant for target
    target_variant = []
    target_supp_reads = partition.supp_reads
    target_supp_num = len(target_supp_reads)

    # # detect from detail cigar
    for detail_cigar in partition.included_hyper_cigars[0].detail_cigars:

        if target_interval_coverage == []:
            detect_pseudo_vaf_for_detail_cigar_via_fetch_bam(detail_cigar, target_supp_num, options.target_path, options.min_mapq)
        else:
            detect_pseudo_vaf_for_detail_cigar_via_coverage_list(detail_cigar, target_supp_num, interval_start, interval_end, target_interval_coverage, options)

        target_variant.append(detect_variant_from_detail_cigar(detail_cigar, target_supp_reads, options))

    # # set info
    partition.set_target_variant(target_variant)


def detect_base_pseudo_variant_for_partition(partition, partition_start, partition_end, base_index, base_interval_coverage, options):
    """
    segment and detect SV from single read's non_match_cigar
    one target variant may have more than one base variants
    """

    # # traverse each base partition
    for base_partition in partition.base_partitions[base_index]:
        cur_base_variant = []
        cur_base_supp_reads = base_partition.supp_reads
        cur_base_supp_num = len(cur_base_supp_reads)

        # # detect from detail cigar
        for detail_cigar in base_partition.included_hyper_cigars[0].detail_cigars:
            if base_interval_coverage == []:
                detect_pseudo_vaf_for_detail_cigar_via_fetch_bam(detail_cigar, cur_base_supp_num, options.base_path, options.min_mapq)
            else:
                detect_pseudo_vaf_for_detail_cigar_via_coverage_list(detail_cigar, cur_base_supp_num, partition_start, partition_end, base_interval_coverage, options)

            cur_base_supp_reads.append(detect_variant_from_detail_cigar(detail_cigar, cur_base_supp_reads, options))

        # # set info
        partition.set_base_variants(base_index, cur_base_variant)

