import re
import numpy as np
import pysam
from src.align_op import fetch_ref_seq, collect_and_boost_supp_aligns, Align, update_coverage_by_aligns, cal_align_mismatch_ratio



class DetailCigar:
    def __init__(self, op, op_len, ref_chrom, ref_start, ref_end, strand):
        self.op = op
        self.hybrid_length = op_len

        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.strand = strand

        self.secondary_ref_chrom = None
        self.secondary_ref_start = None
        self.secondary_ref_end = None
        self.secondary_ref_strand = None

        self.bnd_ref_chrom = None
        self.bnd_ref_start = None
        self.bnd_ref_end = None
        self.bnd_ref_strand = None

        self.alt_bases = ""

        self.uncovered_flag = False
        self.shorter_flag = False

        self.pseudo_vaf = 0
        self.genuine_vaf = 0
        self.coverage = 1e-12

    def set_pseudo_vaf(self, vaf):
        self.pseudo_vaf = vaf
        self.genuine_vaf = self.pseudo_vaf

    def set_coverage(self, coverage):
        self.coverage = coverage

    def set_uncovered_flag(self, flag):
        self.uncovered_flag = flag

    def set_shorter_flag(self, flag):
        self.shorter_flag = flag

    def set_alt_bases(self, alt_bases):
        self.alt_bases = alt_bases

    def set_bnd_mapping(self, chrom, start, end, strand):
        self.bnd_ref_chrom = chrom
        self.bnd_ref_start = start
        self.bnd_ref_end = end
        self.bnd_ref_strand = strand

    def set_secondary_mapping(self, chrom, start, end, strand):
        self.secondary_ref_chrom = chrom
        self.secondary_ref_start = start
        self.secondary_ref_end = end
        self.secondary_ref_strand = strand

    def set_plot_info(self, ref_pointer, read_pointer, resize_ratio):
        self.plot_ref_pointer = ref_pointer
        self.plot_read_pointer = read_pointer
        self.resize_ratio = resize_ratio

    def to_string(self):
        if self.secondary_ref_chrom is not None:
            secondary_info = "{},{},{},{}".format(self.secondary_ref_chrom, self.secondary_ref_start, self.secondary_ref_end, self.secondary_ref_strand)
        else:
            secondary_info = "None"

        if self.bnd_ref_chrom is not None:
            bnd_info = "{},{},{},{}".format(self.bnd_ref_chrom, self.bnd_ref_start, self.bnd_ref_end, self.bnd_ref_strand)
        else:
            bnd_info = "None"

        return "{},{},{},{},{},{},secondary:{},bnd:{}".format(self.op, self.hybrid_length, self.ref_chrom, self.ref_start, self.ref_end, self.strand, secondary_info, bnd_info)
        # return "{},{},{},{},{},{},secondary:{},bnd:{},seq:{}".format(self.op, self.hybrid_length, self.ref_chrom, self.ref_start, self.ref_end, self.strand, secondary_info, bnd_info, self.alt_bases)


class HyperCigar:
    def __init__(self, detail_cigars, supp_reads, uncovered_flag=False, shorter_flag=False):

        self.detail_cigars = detail_cigars

        self.alt_bases = ""
        for cigar in detail_cigars:
            self.alt_bases += cigar.alt_bases

        self.op = "+".join([detail_cigar.op if detail_cigar.op not in ["MappedI", "UnmappedI"] else "I" for detail_cigar in self.detail_cigars])

        self.ref_chrom = self.detail_cigars[0].ref_chrom
        self.ref_start = min([detail_cigar.ref_start for detail_cigar in self.detail_cigars])
        self.ref_end = max([detail_cigar.ref_end for detail_cigar in self.detail_cigars])

        included_chroms = [detail_cigar.ref_chrom for detail_cigar in self.detail_cigars]
        included_chroms.extend([detail_cigar.secondary_ref_chrom for detail_cigar in self.detail_cigars if detail_cigar.secondary_ref_chrom is not None])

        if len(set(included_chroms)) == 1:
            self.ref_start_with_secondary = min([detail_cigar.ref_start if detail_cigar.secondary_ref_start is None else min(detail_cigar.ref_start, detail_cigar.secondary_ref_start) for detail_cigar in self.detail_cigars])
            self.ref_end_with_secondary = max([detail_cigar.ref_end if detail_cigar.secondary_ref_end is None else max(detail_cigar.ref_start, detail_cigar.secondary_ref_end) for detail_cigar in self.detail_cigars])
        else:
            self.ref_start_with_secondary = 0
            self.ref_end_with_secondary = 9999999

        self.hybrid_length = sum([detail_cigar.hybrid_length for detail_cigar in self.detail_cigars])
        self.read_length = self.calculate_read_length()
        self.ref_length = self.ref_end - self.ref_start

        self.supp_reads = supp_reads
        self.uncovered_flag = uncovered_flag
        self.shorter_flag = shorter_flag

    def get_cigar_identifier(self):
        self.identifier = "{}-{}-{}-{}-{}-{}-{}-{}".format(self.ref_chrom, self.ref_start, self.ref_end, self.op, self.hybrid_length, self.ref_start_with_secondary, self.ref_end_with_secondary, self.uncovered_flag)

        return self.identifier

    def calculate_read_length(self):
        read_length = 0

        for detail_cigar in self.detail_cigars:
            # # only D dose not need to increase read length
            if detail_cigar.op != "D":
                read_length += self.hybrid_length

        return read_length

    def update_alt_bases(self):
        self.alt_bases = ""
        for cigar in self.detail_cigars:
            self.alt_bases += cigar.alt_bases

    def calculate_pseudo_vaf(self):
        return round(np.average([detail_cigar.pseudo_vaf for detail_cigar in self.detail_cigars]), 2)

    def calculate_coverage(self):
        return int(np.average([detail_cigar.coverage for detail_cigar in self.detail_cigars]))

    def to_string(self, output_read=False):
        if output_read:
            return "{};{}-{}-{};op={};len={};supp={};details={}".format(self.supp_reads, self.ref_chrom, self.ref_start, self.ref_end, self.op, self.hybrid_length, len(self.supp_reads), [detail.to_string() for detail in self.detail_cigars])
        else:
            return "{}-{}-{};op={};len={};supp={};details={}".format(self.ref_chrom, self.ref_start, self.ref_end, self.op, self.hybrid_length, len(self.supp_reads), [detail.to_string() for detail in self.detail_cigars])


# def CIGAR_copy_constructor(cigar, reset_supp_reads=False):
#     if reset_supp_reads:
#         new_cigar = HyperCigar(cigar.op, cigar.hybrid_length, cigar.ref_bases, cigar.alt_bases, cigar.ref_chrom, cigar.ref_start, cigar.ref_end, [], cigar.source)
#         new_cigar.mapped_Is = cigar.mapped_Is.copy() if cigar.mapped_Is is not None else None
#         return new_cigar
#
#     else:
#         new_cigar = HyperCigar(cigar.op, cigar.hybrid_length, cigar.ref_bases, cigar.alt_bases, cigar.ref_chrom, cigar.ref_start, cigar.ref_end, cigar.supp_reads.copy(), cigar.source)
#         new_cigar.mapped_Is = cigar.mapped_Is.copy() if cigar.mapped_Is is not None else None
#         return new_cigar


def cigar_to_list(cigar, rm_clip=True):
    """
    convert cigar string to list
    :return:
    """

    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    if rm_clip == True:
        if ops[0] == "S" or ops[0] == "H":
            ops = ops[1:]
            lengths = lengths[1:]

        if len(ops) > 1:    # # fix on v1.9 for VACmap aligner, there are aligns with CIGAR like '13520S' and there will be no index -1
            if ops[-1] == "S" or ops[-1] == "H":
                ops = ops[: -1]
                lengths = lengths[: -1]

    return ops, lengths


def collect_cigar_from_ref(ref_path, partial_chrom, partial_start, partial_end, options):
    """
    collect from ref
    :param ref_file:
    :return:
    """

    # # STEP: fetch ref seq
    # # NOTE: -1: convert 1-based to 0based, +1: convert [, ) to [, ]
    ref_seq = pysam.FastaFile(ref_path).fetch(partial_chrom, partial_start - 1, partial_end - 1 + 1)

    # # STEP: collect cigar,
    # # we call this non-match, only for consistent with detecting from bam, since this in ref, so only "=" ops, but we mark this as "X"
    non_match_cigars = []

    for pos in range(partial_start, partial_end + 1):
        ref_base = ref_seq[pos - partial_start]
        detial_cigars = [DetailCigar(ref_base, 1, partial_chrom, pos, pos, "+")]

        cigar = HyperCigar(detial_cigars, ["ref"])

        non_match_cigars.append(cigar)

    return non_match_cigars


def merge_split_Is(hyper_cigar, max_sv_size, max_ins_split=4):
    """
    for a hyper cigar A, if its detail cigars are all Is, then we consider this as a special condition caused by repeats or TE insertion

    there are two types:
    1. intra chrom duplication. but the source bpk is to far from the inserted bkp
    2. inter chrom duplication. for now ,we could only convert the duplication to insertion (removed)
    3. A large insertion, but was cut to several pieces (massive ins) due to repeats or TE or segmental duplication

    """
    if "D" in hyper_cigar.op:
        return hyper_cigar

    if "V" in hyper_cigar.op:
        return hyper_cigar

    intra_chrom_flag = False
    inter_chrom_flag = False
    massive_ins_flag = False

    total_ins_length = 0
    total_ins_seq = ""
    for detail_cigar in hyper_cigar.detail_cigars:
        total_ins_length += detail_cigar.hybrid_length
        total_ins_seq += detail_cigar.alt_bases

        # # deal with type2: inter chrom duplication
        if detail_cigar.op == "MappedI" and detail_cigar.secondary_ref_chrom != hyper_cigar.ref_chrom:
            inter_chrom_flag = True

        if detail_cigar.op == "MappedI" and detail_cigar.secondary_ref_chrom == hyper_cigar.ref_chrom and abs(detail_cigar.secondary_ref_start - hyper_cigar.ref_start) > max_sv_size:
            intra_chrom_flag = True

    # # deal with type3: massive insertion
    if len(hyper_cigar.op.split("+")) > max_ins_split:
        massive_ins_flag = True

    # if inter_chrom_flag is True or massive_ins_flag is True or intra_chrom_flag is True:
    if massive_ins_flag is True or intra_chrom_flag is True:
        detail_cigar = DetailCigar("UnmappedI", total_ins_length, hyper_cigar.ref_chrom, hyper_cigar.ref_start, hyper_cigar.ref_end, "+")
        detail_cigar.set_alt_bases(total_ins_seq)
        return HyperCigar([detail_cigar], hyper_cigar.supp_reads, uncovered_flag=hyper_cigar.uncovered_flag, shorter_flag=hyper_cigar.shorter_flag)
    else:
        return hyper_cigar


def merge_continuous_Is_and_ds(hyper_cigars, dist_merge_continuous):
    """
    for two near hyper cigar Is and Ds, we merge them as one

    """
    # # STEP: find continuous Is and Ds
    continuous_op_indexes = []

    continuous_start_index = 0
    continuous_ref_pos = hyper_cigars[0].ref_end
    continuous_op = hyper_cigars[0].op

    # # STEP: if cur continuous_op is like I+I+I, we consider it as I
    continuous_op_split = continuous_op.split("+")
    if len(set(continuous_op_split)) == 1 and continuous_op_split[0] == "I":
        continuous_op = "I"

    for cur_index in range(1, len(hyper_cigars)):
        cur_op = hyper_cigars[cur_index].op

        # # STEP: if cur op is like I+I+I, we consider it as I
        cur_op_split = cur_op.split("+")
        if len(set(cur_op_split)) == 1 and cur_op_split[0] == "I":
            cur_op = "I"

        # # STEP: locate at a certain distance
        if abs(hyper_cigars[cur_index].ref_start - continuous_ref_pos) > dist_merge_continuous:
            if cur_index - continuous_start_index > 1:
                continuous_op_indexes.append([continuous_op, continuous_start_index, cur_index])

            continuous_start_index = cur_index
            continuous_op = cur_op
            continuous_ref_pos = hyper_cigars[cur_index].ref_end
        else:

            if cur_op != continuous_op:
                if cur_index - continuous_start_index > 1:
                    continuous_op_indexes.append([continuous_op, continuous_start_index, cur_index])

                continuous_op = cur_op
                continuous_start_index = cur_index
                continuous_ref_pos = hyper_cigars[cur_index].ref_end
            else:
                # # meet the last cigar
                if cur_index == len(hyper_cigars) - 1:
                    if cur_index + 1 - continuous_start_index > 1:
                        continuous_op_indexes.append([continuous_op, continuous_start_index, cur_index + 1])

                    continuous_op = cur_op
                    continuous_start_index = cur_index
                    continuous_ref_pos = hyper_cigars[cur_index].ref_end

                continuous_ref_pos = hyper_cigars[cur_index].ref_end
    # # STEP: merge
    deleted_cigar_index = []
    # print(continuous_op_indexes)

    for op, start_index, end_index in continuous_op_indexes:

        # # # when merging ,if there is a uncovered op, then we donot merge them
        # uncover_flag = False
        # for index in range(start_index, end_index):
        #     if hyper_cigars[index].uncovered_flag is True:
        #         uncover_flag = True
        #         break
        #
        # if uncover_flag:
        #     continue

        # # find the longest D or I
        longest_index = start_index
        for index in range(start_index + 1, end_index):
            if hyper_cigars[index].hybrid_length > hyper_cigars[longest_index].hybrid_length:
                longest_index = index

        if op == "D":
            # # merge into longest D
            for index in range(start_index, end_index):

                if index == longest_index:
                    continue
                # print(1111111)
                hyper_cigars[longest_index].hybrid_length += hyper_cigars[index].hybrid_length
                hyper_cigars[longest_index].detail_cigars[0].hybrid_length += hyper_cigars[index].hybrid_length

                # # left Ds, extend start position
                if index < longest_index:
                    hyper_cigars[longest_index].ref_start -= hyper_cigars[index].hybrid_length
                    hyper_cigars[longest_index].detail_cigars[0].ref_start -= hyper_cigars[index].hybrid_length
                else:
                    hyper_cigars[longest_index].ref_end += hyper_cigars[index].hybrid_length
                    hyper_cigars[longest_index].detail_cigars[0].ref_end += hyper_cigars[index].hybrid_length

                deleted_cigar_index.append(index)

        # # merge others to the longest one
        elif op == "I":
            # uncover_flag = [hyper_cigars[start_index].uncovered_flag]

            # # find the detail index of first UnmappedI
            unmappedI_index = -1
            for detail_cigar_index in range(len(hyper_cigars[start_index].detail_cigars)):
                if hyper_cigars[start_index].detail_cigars[detail_cigar_index].op == "UnmappedI":
                    unmappedI_index = detail_cigar_index
                    break

            # print(unmappedI_index)
            # # if find unmappedI, then merge into it
            if unmappedI_index != -1:
                for index in range(start_index + 1, end_index):
                    # uncover_flag.append(hyper_cigars[index].uncovered_flag)

                    hyper_cigars[start_index].hybrid_length += hyper_cigars[index].hybrid_length
                    hyper_cigars[start_index].detail_cigars[unmappedI_index].hybrid_length += hyper_cigars[index].hybrid_length

                    for detail_cigar_index in range(len(hyper_cigars[index].detail_cigars)):
                        hyper_cigars[start_index].detail_cigars[unmappedI_index].alt_bases += hyper_cigars[index].detail_cigars[detail_cigar_index].alt_bases

                    hyper_cigars[start_index].update_alt_bases()

                    # print(len(hyper_cigars[start_index].alt_bases))
                    deleted_cigar_index.append(index)

            # if True in uncover_flag:
            #     hyper_cigars[start_index].uncovered_flag = True

        else:
            continue
            # # # if not found, then create a detail cigar u
            # else:
            #     for index in range(start_index + 1, end_index):
            #         hyper_cigars[start_index].hybrid_length += hyper_cigars[index].hybrid_length
            #
            #         detail_cigar = DetailCigar("UnmappedI", hyper_cigars[index].hybrid_length,  hyper_cigars[index].ref_chrom,  hyper_cigars[index].ref_start, hyper_cigars[index].ref_end, "+")
            #
            #         detail_cigar.set_ref_and_alt_bases("N", "".join())

        """
        
        # # merge others to the first one
        elif op == "I":
            # # find the detail index of first UnmappedI
            unmappedI_index = -1
            for detail_cigar_index in range(len(hyper_cigars[start_index].detail_cigars)):
                if hyper_cigars[start_index].detail_cigars[detail_cigar_index].op == "UnmappedI":
                    unmappedI_index = detail_cigar_index
                    break

            # # if find unmappedI, then merge into it
            if unmappedI_index != -1:
                for index in range(start_index + 1, end_index):

                    hyper_cigars[start_index].hybrid_length += hyper_cigars[index].hybrid_length
                    hyper_cigars[start_index].detail_cigars[unmappedI_index].hybrid_length += hyper_cigars[index].hybrid_length

                    for detail_cigar_index in range(len(hyper_cigars[index].detail_cigars)):
                        hyper_cigars[start_index].detail_cigars[unmappedI_index].alt_bases += hyper_cigars[index].detail_cigars[detail_cigar_index].alt_bases

                    deleted_cigar_index.append(index)

            # # # if not found, then create a detail cigar u
            # else:
            #     for index in range(start_index + 1, end_index):
            #         hyper_cigars[start_index].hybrid_length += hyper_cigars[index].hybrid_length
            #
            #         detail_cigar = DetailCigar("UnmappedI", hyper_cigars[index].hybrid_length,  hyper_cigars[index].ref_chrom,  hyper_cigars[index].ref_start, hyper_cigars[index].ref_end, "+")
            #
            #         detail_cigar.set_ref_and_alt_bases("N", "".join())
        """


    # # delete
    for index in sorted(deleted_cigar_index, reverse=True):
        hyper_cigars.pop(index)


def collect_cigars_from_bam(bam_path, interval_chrom, interval_start, interval_end, mode, options):
    """
    traverse each reads to collect signatures
    :return:
    """
    partial_bam_file = pysam.AlignmentFile(bam_path).fetch(interval_chrom, interval_start, interval_end)

    hyper_cigars_all_reads = {}

    if options.coverage_method == "coverage_list":
        coverage_list = np.zeros(interval_end - interval_start + 1)
    else:
        coverage_list = []

    for align in partial_bam_file:
        hyper_cigars_single_read = []

        # # no cigar, then pass this align
        if align.cigarstring is None:
            continue
        # # unmapped or secondary or low mapping quality, then pass this align
        if align.is_unmapped:
            continue
        if mode == "target" and align.is_secondary:
            continue
        if mode == "target" and align.mapq < options.min_mapq:
            continue

        # # align to a ref that not in genome reference
        align_chr = align.reference_name
        if align_chr not in options.accessible_chroms:
            continue

        if options.mismatch_filter:

            align_mismatch_ratio = cal_align_mismatch_ratio(align.cigarstring, align.reference_start, align.reference_end)

            if align_mismatch_ratio >= options.mismatch_ratio:
                continue

        # #-----------------------Step 1. Collect align's information-----------------------
        # #

        # #-----------------------Step 2. Traverse intra align-----------------------
        # # traverse intra align to find candidate
        hyper_cigars_from_intra, read_aligns_from_intra = collect_from_intra_align(align, mode, options)
        hyper_cigars_single_read.extend(hyper_cigars_from_intra)

        # for cigar in hyper_cigars_from_intra:
        #     print("intra cigar", cigar.uncovered_flag, cigar.shorter_flag, cigar.to_string(output_read=True), len(cigar.alt_bases))

        # #-----------------------Step 3. Traverse inter align-----------------------
        # # whem meet primary align, we process the whole read (inter aligns)
        if not align.is_supplementary:
            hyper_cigars_from_inter, read_aligns_from_inter = collect_from_inter_align(align, mode, options)

            # for cigar in hyper_cigars_from_inter:
            #     print("inter cigar", cigar.uncovered_flag, cigar.shorter_flag, cigar.to_string(output_read=True))

            # if not options.preset == "asm":
            hyper_cigars_single_read.extend(hyper_cigars_from_inter)

            if coverage_list != []:
                update_coverage_by_aligns(align.qname, coverage_list, interval_chrom, interval_start, interval_end, read_aligns_from_inter, options)

        hyper_cigars_single_read = sorted(hyper_cigars_single_read, key=lambda x: x.ref_start)

        # # STEP: for near hyper cigar Is and Ds, we merge them as one
        # if len(hyper_cigars_single_read) != 0 and not options.preset == "asm":
        if len(hyper_cigars_single_read) != 0 and not options.skip_nearby:
            merge_continuous_Is_and_ds(hyper_cigars_single_read, options.dist_merge_continuous)

        # # STEP: add cur read's non match cigars to all reads
        for cigar in hyper_cigars_single_read:
            if cigar.ref_chrom != interval_chrom:
                continue

            if cigar.hybrid_length > options.max_sv_size:
                continue

            cigar = merge_split_Is(cigar, options.max_sv_size)

            if cigar.get_cigar_identifier() not in hyper_cigars_all_reads:
                hyper_cigars_all_reads[cigar.get_cigar_identifier()] = cigar
            else:
                hyper_cigars_all_reads[cigar.get_cigar_identifier()].supp_reads.extend(cigar.supp_reads)

    # # STEP: filter cigar by min support and sort by ref cords

    if not options.skip_coverage_filter and mode == "target":
        mass_coverage_list = coverage_list.copy()
        mass_coverage_list[coverage_list <= options.max_coverage] = 0
        mass_coverage_list[coverage_list > options.max_coverage] = 1

        min_coverage_index = 0
        max_coverage_index = interval_end - interval_start
    else:
        mass_coverage_list = []
        min_coverage_index = 0
        max_coverage_index = interval_end - interval_start

    filtered_hyper_cigars = []

    for cigar_id in hyper_cigars_all_reads:

        cigar = hyper_cigars_all_reads[cigar_id]

        cigar.supp_reads = list(set(cigar.supp_reads))

        if len(cigar.supp_reads) >= 1:

            if not options.skip_coverage_filter and mode == "target":

                cigar_ref_start_index = cigar.ref_start - interval_start
                cigar_read_end_index = cigar.ref_end - interval_start

                if cigar_ref_start_index > max_coverage_index:
                    cigar_ref_start_index = max_coverage_index
                if cigar_ref_start_index < min_coverage_index:
                    cigar_ref_start_index = min_coverage_index
                if cigar_read_end_index > max_coverage_index:
                    cigar_read_end_index = max_coverage_index
                if cigar_read_end_index < min_coverage_index:
                    cigar_read_end_index = min_coverage_index

                # # filter by mass coverage
                if mass_coverage_list[cigar_ref_start_index] == 0 and mass_coverage_list[cigar_read_end_index] == 0:
                    filtered_hyper_cigars.append(cigar)

            else:
                filtered_hyper_cigars.append(cigar)

    filtered_hyper_cigars = sorted(filtered_hyper_cigars, key=lambda x: x.ref_start)

    return filtered_hyper_cigars, coverage_list


def collect_from_intra_align(align, mode, options):
    """
    process intra align, and search candidate cigar ops

    :return:
    """
    # if mode == "target":
    min_sv_size = options.min_sv_size
    # else:
    #     min_sv_size = options.min_sv_size - 10

    # # STEP: extract align's ref info
    align_qname = align.qname
    align_cigar_str = align.cigarstring
    align_seq = align.query

    align_ref_chrom = align.reference_name
    align_ref_start = align.reference_start + 1
    align_ref_end = align.reference_end
    # align_ref_seq = fetch_ref_seq(align_ref_chrom, align_ref_start, align_ref_end, ref_file)

    if align.is_reverse:
        align_obj = Align(align.reference_name, align.reference_start + 1, align.reference_end, align.query_alignment_start + 1, align.query_alignment_end, "-", "{}V".format(align.reference_end - align.reference_start), "Intra")
    else:
        align_obj = Align(align.reference_name, align.reference_start + 1, align.reference_end, align.query_alignment_start + 1, align.query_alignment_end, "+", "{}M".format(align.reference_end - align.reference_start), "Intra")

    # # STEP: traverse cigar to collect non match cigars
    candidate_hyper_cigars = []

    cigar_ops, cigar_ops_length = cigar_to_list(align_cigar_str, rm_clip=True)

    ref_pointer = 0
    read_pointer = 0

    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        op_len = cigar_ops_length[i]

        if op == "N" or op == "S":
            read_pointer += op_len

        elif op == "I":
            # # base included in this op
            alt_bases = align_seq[read_pointer: read_pointer + op_len]

            # # save to list
            if op_len >= min_sv_size:
                cigar_ref_chrom = align_ref_chrom
                cigar_ref_start = ref_pointer + align_ref_start - 1
                cigar_ref_end = ref_pointer + align_ref_start - 1
                detail_cigar = DetailCigar("UnmappedI", op_len, cigar_ref_chrom, cigar_ref_start, cigar_ref_end, "+")

                detail_cigar.set_alt_bases(alt_bases)

                candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname]))

            elif op_len >= min_sv_size - 20:
                cigar_ref_chrom = align_ref_chrom
                cigar_ref_start = ref_pointer + align_ref_start - 1
                cigar_ref_end = ref_pointer + align_ref_start - 1
                detail_cigar = DetailCigar("UnmappedI", op_len, cigar_ref_chrom, cigar_ref_start, cigar_ref_end, "+")

                if mode == "target":
                    detail_cigar.set_shorter_flag(True)
                    candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname], shorter_flag=True))
                else:
                    candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname]))

            # # save to list
            read_pointer += op_len

        elif op == "D":
            # # update base included in this op
            # alt_bases = align_ref_seq[ref_pointer: ref_pointer + op_len]
            alt_bases = ""
            # # save to list
            if op_len >= min_sv_size:
            # if op_len >= options.min_sv_size and "N" not in ref_bases:
                cigar_ref_chrom = align_ref_chrom
                cigar_ref_start = ref_pointer + align_ref_start
                cigar_ref_end = ref_pointer + align_ref_start + op_len - 1

                detail_cigar = DetailCigar(op, op_len, cigar_ref_chrom, cigar_ref_start, cigar_ref_end, "+")
                detail_cigar.set_alt_bases(alt_bases)

                candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname]))
            elif op_len >= min_sv_size - 20:
            # elif op_len >= options.min_sv_size - 10 and "N" not in ref_bases:
                cigar_ref_chrom = align_ref_chrom
                cigar_ref_start = ref_pointer + align_ref_start - 1
                cigar_ref_end = ref_pointer + align_ref_start + op_len - 1
                detail_cigar = DetailCigar(op, op_len, cigar_ref_chrom, cigar_ref_start, cigar_ref_end, "+")

                if mode == "target":
                    detail_cigar.set_shorter_flag(True)
                    candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname], shorter_flag=True))
                else:
                    candidate_hyper_cigars.append(HyperCigar([detail_cigar], [align_qname]))

            ref_pointer += op_len

        elif op in ["M", '=']:
            ref_pointer += op_len
            read_pointer += op_len

        elif op in ["X", "E"]:
            """ commnent SNP processing since we only focus on SVs now
            # # update base included in this op
            ref_base = align_ref_seq[ref_pointer: ref_pointer + op_len]
            alt_bases = align_seq[read_pointer: read_pointer + op_len]

            # # save to list
            for i in range(op_len):
                X_op = alt_bases[i]
                cigar_ref_chrom = align_ref_name
                cigar_ref_start = ref_pointer + align_ref_start
                cigar_ref_end = ref_pointer + align_ref_start
                cigar_details = [[op, cigar_ref_chrom, cigar_ref_start, cigar_ref_end]]

                non_match_cigars.append(CIGAR(X_op, 1, ref_base, alt_bases, cigar_ref_chrom, cigar_ref_start, cigar_ref_end, [align_qname], cigar_details))

            """

            ref_pointer += op_len
            read_pointer += op_len

        else:
            continue

    return candidate_hyper_cigars, [align_obj]


def collect_from_inter_align(primary, mode, options):
    """
    process inter aligns (between primary and supplementary aligns)
    :return:
    """
    read_name = primary.qname

    # # STEP: first, for each primary align, we need collect its supp aligns
    read_aligns, read_align_chroms = collect_and_boost_supp_aligns(primary, mode, options)
    
    # # STEP: filter abnormal aligns
    pm_index = -1
    abnormal_flag = False
    for index in range(len(read_aligns)):
        align = read_aligns[index]

        # # too long ref span
        align_ref_span = align.ref_end - align.ref_start
        if align_ref_span > options.max_sv_size:
            abnormal_flag = True

        align_type = align.source

        if align_type == "PM":
            pm_index = index

    if abnormal_flag is True:
        return [], [read_aligns[pm_index]]

    # # set seq
    read_seq = primary.query_sequence
    for align in read_aligns:
        # if "D" in align.cigar:
        #     align.set_align_seq(fetch_ref_seq(align.ref_chrom, align.ref_start, align.ref_end, ref_file))
        # else:
            align.set_align_seq(read_seq[align.read_start: align.read_end])

    # # STEP: traverse between match aligns to find non match aligns
    candidate_hyper_cigars = []

    # # STEP: inter-chrom, BNDs
    if not options.skip_bnd:
        for align_index in range(len(read_aligns) - 1):
            cur_align = read_aligns[align_index]
            next_align = read_aligns[align_index + 1]

            # # different chrom, then we create a BND cigar
            if cur_align.ref_chrom != next_align.ref_chrom:

                # if cur_align.ref_chrom < next_align.ref_chrom:

                if cur_align.strand == "+":
                    tmp_detial_cigar = DetailCigar("B", options.min_sv_size, cur_align.ref_chrom, cur_align.ref_end, cur_align.ref_end, cur_align.strand)
                else:
                    tmp_detial_cigar = DetailCigar("B", options.min_sv_size, cur_align.ref_chrom, cur_align.ref_start, cur_align.ref_start, cur_align.strand)

                if next_align.strand == "+":
                    tmp_detial_cigar.set_bnd_mapping(next_align.ref_chrom, next_align.ref_start, next_align.ref_start, next_align.strand)
                else:
                    tmp_detial_cigar.set_bnd_mapping(next_align.ref_chrom, next_align.ref_end, next_align.ref_end, next_align.strand)

                # else:
                #     tmp_detial_cigar = DetailCigar("B", -1, next_align.ref_chrom, next_align.ref_start, next_align.ref_start, next_align.strand)
                #     tmp_detial_cigar.set_bnd_mapping(cur_align.ref_chrom, cur_align.ref_end, cur_align.ref_end, cur_align.strand)

                candidate_hyper_cigars.append(HyperCigar([tmp_detial_cigar], [read_name]))

    # # STEP: intra-chrom, BND and non-BND events
    for spec_chrom in read_align_chroms:
        read_aligns_spec_chrom = [align for align in read_aligns if align.ref_chrom == spec_chrom]

        # # STEP: find match aligns, store their index
        match_aligns_index = [i for i in range(len(read_aligns_spec_chrom)) if "M" in read_aligns_spec_chrom[i].cigar]

        # # STEP: aligns between each two match aligns are non_match aligns
        for i in range(len(match_aligns_index) - 1):
            cur_match_align_index = match_aligns_index[i]
            next_match_align_index = match_aligns_index[i + 1]

            included_non_match_aligns = read_aligns_spec_chrom[cur_match_align_index + 1: next_match_align_index]
            included_non_match_aligns_chroms = list(set([ali.ref_chrom for ali in included_non_match_aligns]))

            if len(included_non_match_aligns) == 0:
                continue
            if len(included_non_match_aligns_chroms) != 1:
                continue
            if included_non_match_aligns_chroms[0] != spec_chrom:
                continue

            # # STEP: generate hyper cigars.
            hyper_cigar_ops = []
            hyper_cigar_op_lens = []
            hyper_cigar_details = []

            # # STEP: each align represents a detail cigar.
            for non_match_align in included_non_match_aligns:
                align_cigar_op, align_cigar_op_len = cigar_to_list(non_match_align.cigar)

                if len(align_cigar_op) == 0:     # # fix on v1.9 for empty align_cigar_op
                    continue

                if options.skip_flanking_indels and align_cigar_op_len[0] < options.min_sv_size:
                    continue

                # # Is are special, since they could represent duplication, therefore, we mark unmapped or mapped Is
                if align_cigar_op[0] == "I":
                    if non_match_align.secondary_ref_chrom is None:
                        tmp_detial_cigar = DetailCigar("UnmappedI", align_cigar_op_len[0], non_match_align.ref_chrom, non_match_align.ref_start, non_match_align.ref_end, non_match_align.strand)
                        tmp_detial_cigar.set_alt_bases(non_match_align.seq)

                    else:
                        # # # different chroms mean inter-chrom translocation, for now, we do not detect inter events
                        # if non_match_align.secondary_ref_chrom != non_match_align.ref_chrom:
                        #     tmp_detial_cigar = DetailCigar("UnmappedI", align_cigar_op_len[0], non_match_align.ref_chrom, non_match_align.ref_start, non_match_align.ref_end, non_match_align.strand)
                        #     tmp_detial_cigar.set_alt_bases(non_match_align.seq)
                        #
                        # else:
                        # # set mapping info
                        tmp_detial_cigar = DetailCigar("MappedI", align_cigar_op_len[0], non_match_align.ref_chrom, non_match_align.ref_start, non_match_align.ref_end, non_match_align.strand)
                        tmp_detial_cigar.set_alt_bases(non_match_align.seq)
                        tmp_detial_cigar.set_secondary_mapping(non_match_align.secondary_ref_chrom, non_match_align.secondary_ref_start, non_match_align.secondary_ref_end, non_match_align.secondary_strand)

                else:
                    # # intra-chrom BNDs
                    if non_match_align.ref_end - non_match_align.ref_start > options.max_sv_size:
                        tmp_detial_cigar = DetailCigar("B", 1, non_match_align.ref_chrom, non_match_align.ref_start, non_match_align.ref_start, non_match_align.strand)
                        tmp_detial_cigar.set_bnd_mapping(non_match_align.ref_chrom, non_match_align.ref_end, non_match_align.ref_end, non_match_align.strand)
                    # # for D and V
                    else:
                        # if "N" in non_match_align.seq:
                        #     continue
                        tmp_detial_cigar = DetailCigar(align_cigar_op[0], align_cigar_op_len[0], non_match_align.ref_chrom, non_match_align.ref_start, non_match_align.ref_end, non_match_align.strand)
                        tmp_detial_cigar.set_alt_bases(non_match_align.seq)

                # # append
                if options.skip_bnd and tmp_detial_cigar.op == "B":
                    continue

                hyper_cigar_ops.append(align_cigar_op[0])
                hyper_cigar_op_lens.append(align_cigar_op_len[0])
                hyper_cigar_details.append(tmp_detial_cigar)

            if len(hyper_cigar_ops) != 0:
                candidate_hyper_cigars.append(HyperCigar(hyper_cigar_details, [read_name]))

    # if options.detect_mode == "germline":
    # # STEP: check whether the left most aligns are I
    hyper_cigar_details = []
    for i in range(len(read_aligns)):
        tmp_align = read_aligns[i]
        tmp_align_cigar_op, tmp_align_cigar_op_len = cigar_to_list(tmp_align.cigar)

        if len(tmp_align_cigar_op) == 0:    # # fix on v1.9 for empty align_cigar_op
            continue

        if tmp_align_cigar_op[0] != "I":
            break
        else:
            if tmp_align.secondary_ref_chrom is not None:
                tmp_detial_cigar = DetailCigar("MappedI", tmp_align_cigar_op_len[0], tmp_align.ref_chrom, tmp_align.ref_start, tmp_align.ref_end, tmp_align.strand)
                tmp_detial_cigar.set_secondary_mapping(tmp_align.secondary_ref_chrom, tmp_align.secondary_ref_start, tmp_align.secondary_ref_end, tmp_align.secondary_strand)
            else:
                tmp_detial_cigar = DetailCigar("UnmappedI", tmp_align_cigar_op_len[0], tmp_align.ref_chrom, tmp_align.ref_start, tmp_align.ref_end, tmp_align.strand)
            tmp_detial_cigar.set_uncovered_flag(True)
            tmp_detial_cigar.set_alt_bases(tmp_align.seq)
            hyper_cigar_details.append(tmp_detial_cigar)

    if len(hyper_cigar_details) != 0:
        candidate_hyper_cigars.append(HyperCigar(hyper_cigar_details, [read_name], uncovered_flag=True))

    # # STEP: check whether the right most aligns are I
    hyper_cigar_details = []
    for i in range(len(read_aligns) - 1, -1, -1):
        tmp_align = read_aligns[i]
        tmp_align_cigar_op, tmp_align_cigar_op_len = cigar_to_list(tmp_align.cigar)

        if len(tmp_align_cigar_op) == 0:    # # fix on v1.9 for empty align_cigar_op
            continue

        if tmp_align_cigar_op[0] != "I":
            break
        else:
            if tmp_align.secondary_ref_chrom is not None:
                tmp_detial_cigar = DetailCigar("MappedI", tmp_align_cigar_op_len[0], tmp_align.ref_chrom, tmp_align.ref_start, tmp_align.ref_end, tmp_align.strand)
                tmp_detial_cigar.set_secondary_mapping(tmp_align.secondary_ref_chrom, tmp_align.secondary_ref_start, tmp_align.secondary_ref_end, tmp_align.secondary_strand)
            else:
                tmp_detial_cigar = DetailCigar("UnmappedI", tmp_align_cigar_op_len[0], tmp_align.ref_chrom, tmp_align.ref_start, tmp_align.ref_end, tmp_align.strand)

            # tmp_detial_cigar.set_uncovered_flag(True)
            tmp_detial_cigar.set_alt_bases(tmp_align.seq)

            hyper_cigar_details.append(tmp_detial_cigar)

    if len(hyper_cigar_details) != 0:
        candidate_hyper_cigars.append(HyperCigar(hyper_cigar_details, [read_name], uncovered_flag=True))

    # for cigar in candidate_hyper_cigars:
    #     print(1111, cigar.uncovered_flag, cigar.to_string())
    return candidate_hyper_cigars, read_aligns
