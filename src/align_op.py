import re
from collections import Counter


class Align:
    def __init__(self, ref_chrom, ref_start, ref_end, read_start, read_end, strand, cigar, source):
        self.ref_chrom = ref_chrom
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end

        self.origin_ref_chrom = ref_chrom
        self.origin_ref_start = ref_start
        self.origin_ref_end = ref_end

        self.strand = strand
        self.cigar = cigar

        self.source = source
        self.secondary_ref_chrom = None
        self.secondary_ref_start = None
        self.secondary_ref_end = None
        self.secondary_strand = None
        self.seq = ""
        self.id = None

    def set_align_id(self, id):
        self.id = id

    def set_secondary_mapping(self, ref_chrom, ref_start, ref_end, strand):
        self.secondary_ref_chrom = ref_chrom
        self.secondary_ref_start = ref_start
        self.secondary_ref_end = ref_end
        self.secondary_strand = strand

    def set_align_seq(self, seq):
        self.seq = seq

    def to_string(self):
        if self.secondary_ref_chrom is not None:
            return "{}-{}-{}-{}-{}-{}-{}-{}-{}-{}-{}".format(self.source, self.ref_chrom, self.ref_start, self.ref_end, self.read_start, self.read_end, self.strand, self.cigar, self.secondary_ref_chrom, self.secondary_ref_start, self.secondary_ref_end)
        else:
            return "{}-{}-{}-{}-{}-{}-{}-{}".format(self.source, self.ref_chrom, self.ref_start, self.ref_end, self.read_start, self.read_end, self.strand, self.cigar)


def calculate_ref_and_read_end_by_cigar(ref_start, read_start, cigar_ops, cigar_ops_length):
    """
    given ref and read start, calculate ref and read end by cigar ops
    :param ref_start:
    :param read_start:
    :param ops:
    :param ops_len:
    :return:
    """

    ref_end = ref_start
    read_end = read_start
    for i in range(len(cigar_ops)):
        op = cigar_ops[i]
        op_len = cigar_ops_length[i]

        if op == "N":
            read_end += op_len

        elif op == "I":
            read_end += op_len

        elif op == "D":
            ref_end += op_len

        elif op in ["M", '=', "X", "E"]:
            ref_end += op_len
            read_end += op_len
        else:
            continue

    return ref_end - 1, read_end - 1


def fetch_ref_seq(chrom, start, end, ref_file):
    """
    fetch ref sequence by the chrom, start and end
    :return:
    """
    # -1: convert 1-based to 0-based; +1: convert [ , ) to [ , ]
    return ref_file.fetch(chrom, start - 1, end - 1 + 1)


def update_coverage_by_aligns(read_name, coverage_list, list_ref_chrom, list_ref_start, list_ref_end, align_obj_list, options):
    """
    updata coverage list by aligns, add the covered position by one
    """

    # # STEP: collect start end position pairs for each align
    ref_start_end_pairs = []

    for align_obj in align_obj_list:
        # print(align_obj.to_string())
        # # not same chrom
        if align_obj.ref_chrom != list_ref_chrom:
            continue

        # # only care Intra and D aligns, since we do not want count M align twice
        # if align_obj.source not in ["Intra", "D"]:
        #     continue

        # # skip I align, only D, V, M aligns
        if "I" in align_obj.cigar:
            continue

        if "D" in align_obj.cigar and align_obj.origin_ref_end - align_obj.origin_ref_start > options.max_sv_size:
            continue

        if not (align_obj.origin_ref_start <= list_ref_start <= align_obj.origin_ref_end or
                align_obj.origin_ref_start <= list_ref_end <= align_obj.origin_ref_end or
                list_ref_start <= align_obj.origin_ref_start <= list_ref_end or
                list_ref_start <= align_obj.origin_ref_end <= list_ref_end):
            continue
        #if not list_ref_start <= align_obj.origin_ref_start <= align_obj.origin_ref_end <= list_ref_end:
        #    continue

        if [align_obj.origin_ref_start, align_obj.origin_ref_end] not in ref_start_end_pairs:
            ref_start_end_pairs.append([align_obj.origin_ref_start, align_obj.origin_ref_end])

    # # STEP: resolve overlapped aligns and only increase once (unique read)
    ref_start_end_pairs = sorted(ref_start_end_pairs, key=lambda x: x[0])
    ref_start_end_pairs_final = []

    i = 0
    while i <= len(ref_start_end_pairs) - 1:
        if i == len(ref_start_end_pairs) - 1:
            ref_start_end_pairs_final.append([ref_start_end_pairs[i][0], ref_start_end_pairs[i][1]])
        else:
            cur_ref_start, cur_ref_end = ref_start_end_pairs[i][0], ref_start_end_pairs[i][1]
            next_ref_start, next_ref_end = ref_start_end_pairs[i + 1][0], ref_start_end_pairs[i + 1][1]

            ref_start_end_pairs_final.append([cur_ref_start, cur_ref_end])

            if next_ref_start <= cur_ref_end:
                # # fully overlapped, skip this align
                if next_ref_end <= cur_ref_end:
                    i += 1
                # # otherwise, update this align's start pos
                else:
                    ref_start_end_pairs[i + 1][0] = cur_ref_end + 1

        i += 1

    for ref_start, ref_end in ref_start_end_pairs_final:
        # coverage_list[ref_start - list_ref_start: ref_end - list_ref_start + 1] += 1
        coverage_list[max(list_ref_start, ref_start) - list_ref_start: min(list_ref_end, ref_end) - list_ref_start + 1] += 1


def collect_and_boost_supp_aligns(primary, mode, options):
    """
    for each primary align, collect all supplementary aligns and boost them (resolve dup and add gaps)
    :return:
    """
    primary_forward = "-" if primary.is_reverse else "+"

    # # STEP: get info of supplementary aligns
    supp_aligns = get_supplementary_aligns(primary, primary_forward, mode, options)

    # if len(supp_aligns) > options.max_split_num:
    #     return [], []

    # # STEP: generate info of primary align
    pm_length = primary.reference_end - primary.reference_start + 1
    pm_align = Align(primary.reference_name, primary.reference_start + 1, primary.reference_end, primary.query_alignment_start + 1, primary.query_alignment_end, "+", "{}M".format(pm_length), "PM")

    # # STEP: collect signatures from pm and sa
    read_aligns = sorted([pm_align] + supp_aligns, key=lambda x: (x.read_start + x.read_end) / 2)

    read_length = primary.query_length

    pm_align_index = -1
    sa_align_num_same_pm_chrom = 0
    for i in range(len(read_aligns)):
        if read_aligns[i].source == "PM":
            pm_align_index = i
        else:
            if read_aligns[i].ref_chrom == pm_align.ref_chrom:
                sa_align_num_same_pm_chrom += 1

    # # FIX, =====================================================================================
    # # when the pm align located at the right most or left most index, and the other aligns are all from other chroms
    if pm_align_index == len(read_aligns) - 1 and sa_align_num_same_pm_chrom == 0:
        for align in read_aligns[0: pm_align_index]:
            align_cigar_op, align_cigar_op_len = cigar_to_list(align.cigar)
            align_cigar_op[0] = "I"
            align.cigar = "{}{}".format(align_cigar_op_len[0], align_cigar_op[0])
            # # update secondary mapping
            align.set_secondary_mapping(align.ref_chrom, align.ref_start, align.ref_end, align.strand)
            # # update ref info to last align
            align.ref_chrom = read_aligns[pm_align_index].ref_chrom
            align.ref_start = read_aligns[pm_align_index].ref_start
            align.ref_end = read_aligns[pm_align_index].ref_start
            align.strand = read_aligns[pm_align_index].strand

    if pm_align_index == 0 and sa_align_num_same_pm_chrom == 0:
        for align in read_aligns[1: ]:
            align_cigar_op, align_cigar_op_len = cigar_to_list(align.cigar)
            align_cigar_op[0] = "I"
            align.cigar = "{}{}".format(align_cigar_op_len[0], align_cigar_op[0])
            # # update secondary mapping
            align.set_secondary_mapping(align.ref_chrom, align.ref_start, align.ref_end, align.strand)
            # # update ref info to last align
            align.ref_chrom = read_aligns[pm_align_index].ref_chrom
            align.ref_start = read_aligns[pm_align_index].ref_end
            align.ref_end = read_aligns[pm_align_index].ref_end
            align.strand = read_aligns[pm_align_index].strand
    # # END. =====================================================================================

    # # STEP: specific primary chrom as the loftmost align's chrom, which was used to generate inserted aligns
    # # this is tricky, since the primary align is longest align to be determined as primary but might be inserted align (such as the insertion of a long TE)
    leftmost_align = read_aligns[0]
    appointed_primary_chrom = leftmost_align.ref_chrom
    if appointed_primary_chrom != primary.reference_name:
        # # adjust all aligns when left most align is in different chrom and reversed
        if leftmost_align.strand == "-":
            for align in read_aligns:
                new_read_start = max(0, read_length - align.read_end)
                new_read_end = min(read_length, read_length - align.read_start)

                align.read_start = new_read_start
                align.read_end = new_read_end

                align_cigar_op, align_cigar_op_len = cigar_to_list(align.cigar)

                if align.strand == "-":
                    align.strand = "+"
                    align_cigar_op[0] = "M"
                    align.cigar = "{}{}".format(align_cigar_op_len[0], align_cigar_op[0])

                else:
                    align.strand = "-"
                    align_cigar_op[0] = "V"
                    align.cigar = "{}{}".format(align_cigar_op_len[0], align_cigar_op[0])
        else:
            pass

    read_aligns = sorted(read_aligns, key=lambda x: (x.read_start + x.read_end) / 2)

    read_align_chroms = []
    for i in range(len(read_aligns)):
        read_aligns[i].set_align_id(i)
        if read_aligns[i].ref_chrom not in read_align_chroms:
            read_align_chroms.append(read_aligns[i].ref_chrom)

    # # STEP: find duplicated aligns
    read_aligns = resolve_duplicated_aligns(read_aligns, read_align_chroms)

    # # STEP: determine BNDs, there are two kinds of BND, single-end BND or double-ends BND
    # # NOTICE: BND are defined as chr_small to chr_large when inter, chr_small_pos to chr_large_pos when intra
    # # this aims to change two-ends BND into dDUPs
    if len(read_align_chroms) > 1:
        # # init left anchor to the first align
        left_anchor_align_index = 0

        for right_anchor_align_index in range(1, len(read_aligns)):

            if read_aligns[right_anchor_align_index].ref_chrom == appointed_primary_chrom:
                # # there is inter-chrom aligns between the two anchors
                if right_anchor_align_index - left_anchor_align_index > 1:
                    # # then we need to change aligns, which are between left and right anchor aligns, to I aligns
                    for align in read_aligns[left_anchor_align_index + 1: right_anchor_align_index]:
                        align_cigar_op, align_cigar_op_len = cigar_to_list(align.cigar)
                        align_cigar_op[0] = "I"
                        align.cigar = "{}{}".format(align_cigar_op_len[0], align_cigar_op[0])
                        # # update secondary mapping
                        align.set_secondary_mapping(align.ref_chrom, align.ref_start, align.ref_end, align.strand)
                        # # update ref info to last align
                        align.ref_chrom = read_aligns[left_anchor_align_index].ref_chrom
                        align.ref_start = read_aligns[left_anchor_align_index].ref_end
                        align.ref_end = read_aligns[left_anchor_align_index].ref_end
                        align.strand = read_aligns[left_anchor_align_index].strand
                left_anchor_align_index = right_anchor_align_index

    # # # STEP: boost read aligns to find unmapped sequence (as novel "I") and deleted sequence (as "D")
    # # STEP: find unmapped sequence
    # # condition 1: unmapped at the left most
    leftmost_align = read_aligns[0]
    insert_length = leftmost_align.read_start - 0
    if options.min_sv_size < insert_length:
        read_aligns.insert(0, Align(primary.reference_name, leftmost_align.ref_start, leftmost_align.ref_start, 0,  leftmost_align.read_start - 1, "+", "{}I".format(insert_length), "Boost"))

    # # condition 2: unmapped at the right most
    rightmost_align = read_aligns[-1]
    insert_length = read_length - rightmost_align.read_end
    if options.min_sv_size < insert_length:
        read_aligns.append(Align(primary.reference_name, rightmost_align.ref_end, rightmost_align.ref_end, rightmost_align.read_end + 1, read_length, "+", "{}I".format(insert_length), "Boost"))

    # # condition 3: unmapped between aligns
    for i in range(len(read_aligns) - 1):
        cur_align_read_end = read_aligns[i].read_end + 1
        next_align_read_start = read_aligns[i + 1].read_start - 1
        insert_length = next_align_read_start - cur_align_read_end + 1
        if options.min_sv_size < insert_length:
            read_aligns.append(Align(primary.reference_name, read_aligns[i].ref_end, read_aligns[i].ref_end, cur_align_read_end, next_align_read_start, "+", "{}I".format(insert_length), "Boost"))

    read_aligns = sorted(read_aligns, key=lambda x: (x.read_start + x.read_end) / 2)

    # STEP: find deleted sequence between same chrom aligns
    for chrom in read_align_chroms:
        read_aligns_spec_chrom = [align for align in read_aligns if align.ref_chrom == chrom]
        for i in range(len(read_aligns_spec_chrom) - 1):
            cur_align_ref_end = read_aligns_spec_chrom[i].ref_end + 1
            next_align_ref_start = read_aligns_spec_chrom[i + 1].ref_start - 1
            deleted_length = next_align_ref_start - cur_align_ref_end + 1

            if options.min_sv_size < deleted_length:
                read_aligns.append(Align(read_aligns_spec_chrom[i].ref_chrom, cur_align_ref_end, next_align_ref_start, read_aligns_spec_chrom[i].read_end, read_aligns_spec_chrom[i].read_end, "+", "{}D".format(deleted_length), "Boost"))

    read_aligns = sorted(read_aligns, key=lambda x: (x.read_start + x.read_end) / 2)

    # resolve_linear_ins_del(read_aligns, options)

    return read_aligns, read_align_chroms


def get_supplementary_aligns(align, pm_strand, mode, options):
    """
    for primary aligns, obtain their supplementary aligns through SA tag

    """

    try:
        sa_tags = align.get_tag("SA").split(";")
    except KeyError:
        return []

    supplementary_aligns = []

    for tag in sa_tags:

        tag_split = tag.split(",")
        align_chrom = tag_split[0]

        if align_chrom not in options.accessible_chroms:
            continue

        # print(tag)
        # STEP: remove incomplete info fields
        if len(tag_split) != 6:
            continue

        # print(tag)
        mapq = int(tag_split[4])

        # # STEP: remove when mapq less than min_mapq
        if mode == "target" and mapq < options.min_mapq:
            continue

        # # STEP: calculate align's read start
        cigar_str = tag_split[3]
        align_strand = tag_split[2]

        # obtain soft clips' length
        ops, ops_len = cigar_to_list(cigar_str, rm_clip=False)
        left_soft_clip_length = ops_len[0] if ops[0] == "S" else 0
        right_soft_clip_length = ops_len[-1] if ops[-1] == "S" else 0

        if align_strand == pm_strand:
            align_read_start = left_soft_clip_length + 1
            align_strand = "+"
        else:
            align_read_start = right_soft_clip_length + 1
            align_strand = "-"

        # # STEP: calculate align's reference end
        align_ref_start = int(tag_split[1])
        align_ref_end, align_read_end = calculate_ref_and_read_end_by_cigar(align_ref_start, align_read_start, ops, ops_len)

        sa_mismatch_ratio = int(tag_split[-1]) / (align_ref_end - align_ref_start)
        if options.mismatch_filter and sa_mismatch_ratio >= options.mismatch_ratio:
        # if options.mismatch_filter and sa_mismatch_ratio >= options.mismatch_ratio and options.detect_mode != "germline":
            continue
        # # STEP: generate sa info
        align_length = align_ref_end - align_ref_start + 1
        align_simplified_cigar = "{}M".format(align_length) if align_strand == "+" else "{}V".format(align_length)

        supplementary_aligns.append(Align(align_chrom, align_ref_start, align_ref_end, align_read_start, align_read_end, align_strand, align_simplified_cigar, "SA"))

    return supplementary_aligns


def cal_align_mismatch_ratio(cigar_string, align_ref_start, align_ref_end):
    op_counter = Counter(cigar_string)

    if "X" not in op_counter:
        op_counter["X"] = 0

    if "D" not in op_counter:
        op_counter["D"] = 0

    if "I" not in op_counter:
        op_counter["I"] = 0

    return (op_counter["X"] + op_counter["I"] + op_counter["D"]) / (align_ref_end - align_ref_start)
    # return (op_counter["X"]) / (align_ref_end - align_ref_start)


def cal_overlap_length(A_align, B_align):
    """
    calculate overlap between

    :return:
        whole or parital overlap, and overlap start, end cords
    """

    cord_bias = 50

    if A_align.ref_chrom != B_align.ref_chrom:
        return "non_overlap", -1

    if A_align == B_align:
        return "non_overlap", -1

    # # base is totaly covered by target
    if B_align.ref_start - cord_bias <= A_align.ref_start < A_align.ref_end <= B_align.ref_end + cord_bias:
        overlap_length = A_align.ref_end - A_align.ref_start + 1
        return "A_is_fully_overlapped_by_B", overlap_length

    # # base is totaly covered by target
    if A_align.ref_start - cord_bias <= B_align.ref_start < B_align.ref_end <= A_align.ref_end + cord_bias:
        overlap_length = B_align.ref_end - B_align.ref_start + 1
        return "B_is_fully_overlapped_by_A", overlap_length

    # # base is partially covered by target

    if A_align.ref_start < B_align.ref_start < A_align.ref_end < B_align.ref_end:
        overlap_length = A_align.ref_end - B_align.ref_start + 1
        return "B_is_partially_overlapped_by_A", overlap_length

    if B_align.ref_start < A_align.ref_start < B_align.ref_end < A_align.ref_end:
        overlap_length = B_align.ref_end - A_align.ref_start + 1
        return "A_is_partially_overlapped_by_B", overlap_length

    return "non_overlap", -1


def resolve_linear_ins_del(read_aligns, options):
    """
    for a ins align, when its next (or previous) align is del, and their size are similar, then we consider they are linear and remove them
    """

    align_index = 0

    removed_align_index = []

    while align_index < len(read_aligns) - 1:

        cur_align = read_aligns[align_index]
        next_align = read_aligns[align_index + 1]
        
        cur_align_cigar_op, cur_align_cigar_op_len = cigar_to_list(cur_align.cigar)
        next_align_cigar_op, next_align_cigar_op_len = cigar_to_list(next_align.cigar)

        # # meet close I and D
        if (cur_align_cigar_op[0] == "I" and next_align_cigar_op[0] == "D" and cur_align.secondary_ref_chrom is None) or (cur_align_cigar_op[0] == "D" and next_align_cigar_op[0] == "I" and next_align.secondary_ref_chrom is None):
            # # # directly remove
            # read_aligns.remove(cur_align)
            # read_aligns.remove(next_align)

            # # similar size, remove both
            if abs(cur_align_cigar_op_len[0] - next_align_cigar_op_len[0]) < options.min_sv_size:
                removed_align_index.append(align_index)
                removed_align_index.append(align_index + 1)
            # # otherwise, remove the short one
            else:
                pass
                """
                the following codes, which remove events like INS+DEL (useful for old ONTs that contains many errors) , are only for ONT in early years
                """
                if options.preset == "error-prone":
                    if cur_align_cigar_op_len[0] > next_align_cigar_op_len[0]:
                        # # remove the short
                        removed_align_index.append(align_index + 1)

                        # # update the long
                        subtract_op_length = cur_align_cigar_op_len[0] - next_align_cigar_op_len[0]
                        subtract_op = cur_align_cigar_op[0]

                        cur_align.cigar = "{}{}".format(subtract_op_length, subtract_op)

                        if subtract_op == "D":
                            cur_align.ref_end = cur_align.ref_start + subtract_op_length

                    else:
                        # # remove the short
                        removed_align_index.append(align_index)

                        # # update the long
                        subtract_op_length = next_align_cigar_op_len[0] - cur_align_cigar_op_len[0]
                        subtract_op = next_align_cigar_op[0]

                        next_align.cigar = "{}{}".format(subtract_op_length, subtract_op)

                        if subtract_op == "D":
                            next_align.ref_end = next_align.ref_start + subtract_op_length

            align_index += 1

        align_index += 1

    for index in sorted(removed_align_index, reverse=True):
        read_aligns.pop(index)


def resolve_duplicated_aligns(read_aligns, read_align_chroms):
    """
    resolve overlapped aligns, cut them and generate I
    :param read_aligns:
    :return:
    """

    resolved_read_aligns = read_aligns.copy()

    for chrom in read_align_chroms:
        read_aligns_spec_chrom = [align for align in read_aligns if align.ref_chrom == chrom]

        # # STEP: first, fix the cords of first and last align, make sure they extends to the leftmost and rightmost
        left_most_cord = min([align.ref_start for align in read_aligns_spec_chrom])
        right_most_cord = max([align.ref_end for align in read_aligns_spec_chrom])

        first_align = read_aligns_spec_chrom[0]
        last_align = read_aligns_spec_chrom[-1]

        first_align.ref_start = left_most_cord
        last_align.ref_end = right_most_cord

        # # STEP: set major align as the first align
        major_align = read_aligns_spec_chrom[0]

        # # STEP: traverse the following aligns
        for next_align in read_aligns_spec_chrom[1: ]:

            minor_align_flag = 0

            # # for each base align, find whether is fully covered by others
            for target_align in read_aligns_spec_chrom:
                overlap_type, overlap_length = cal_overlap_length(next_align, target_align)

                # # next_align is fully covered, and therefore there is an insertion
                if overlap_type == "A_is_fully_overlapped_by_B":
                    # print("fully,", "major(A): ", target_align.to_string(), "next: ", next_align.to_string())

                    minor_align_flag = 1
                    insert_pos = major_align.ref_end

                    # # generate new align and append it to the list
                    new_align = Align(major_align.ref_chrom, insert_pos, insert_pos, next_align.read_start, next_align.read_start + overlap_length - 1, "+", "{}I".format(overlap_length), "DUP")
                    new_align.set_secondary_mapping(next_align.ref_chrom, next_align.ref_start, next_align.ref_start + overlap_length - 1, next_align.strand)
                    resolved_read_aligns.append(new_align)
                    resolved_read_aligns.remove(next_align)

                    break

            # # minor flag is not 1, then next align is not fully overlapped by
            if minor_align_flag != 1:

                # # check the overlap condition with major align
                overlap_type, overlap_length = cal_overlap_length(next_align, major_align)

                # # next align is partially overlapped by major align, then there is an insertion
                if overlap_type == "A_is_partially_overlapped_by_B":
                    insert_pos = major_align.ref_end

                    # # generate new align and append it to the list
                    new_align = Align(major_align.ref_chrom, insert_pos, insert_pos, next_align.read_start, next_align.read_start + overlap_length - 1, "+", "{}I".format(overlap_length), "DUP")
                    new_align.set_secondary_mapping(next_align.ref_chrom, next_align.ref_start, next_align.ref_start + overlap_length - 1, next_align.strand)
                    resolved_read_aligns.append(new_align)

                    # # update next align's cords in major_align
                    next_align.ref_start += overlap_length - 1
                    next_align.read_start += overlap_length - 1

                    # print("partially,", "major: ", major_align.to_string(), "next: ", next_align.to_string())

                    # # update major align
                    major_align = next_align

                elif overlap_type == "non_overlap":
                    # # update major align
                    major_align = next_align

                    # print("not overlapped by major, change this to major align")
                    # print("Non overlap,", "major: ", major_align.to_string(), "next: ", next_align.to_string())

                elif overlap_type == "B_is_fully_overlapped_by_A":
                    # # this is a special condition when meeting reads that are not long enough to cover the whole events
                    # # for example a DUP, where the leftmost align is lost since the short read length
                    # # in that case, the major align (B) is fully coverred by next align (A)
                    # # therefore, the insert pos is major_align's ref_start, and recalculate the overlap length to covered next align length
                    insert_pos = major_align.ref_end

                    # # re-calculate overlap length
                    overlap_length = major_align.ref_end - next_align.ref_start

                    # # generate new align and append it to the list
                    new_align = Align(major_align.ref_chrom, insert_pos, insert_pos, next_align.read_start, next_align.read_start + overlap_length - 1, "+", "{}I".format(overlap_length), "DUP")
                    new_align.set_secondary_mapping(next_align.ref_chrom, next_align.ref_start, next_align.ref_start + overlap_length - 1, next_align.strand)
                    resolved_read_aligns.append(new_align)

                    # # update next align's cords in resolved_read_aligns
                    next_align.ref_start += overlap_length - 1
                    next_align.read_start += overlap_length - 1

                    # print("partially,", "major: ", major_align.to_string(), "next: ", next_align.to_string())

                    # # update major align
                    major_align = next_align

                else:
                    insert_pos = major_align.ref_end

                    # # generate new align and append it to the list
                    new_align = Align(major_align.ref_chrom, insert_pos, insert_pos, next_align.read_start, next_align.read_start + overlap_length - 1, "+", "{}I".format(overlap_length), "DUP")
                    new_align.set_secondary_mapping(next_align.ref_chrom, next_align.ref_start, next_align.ref_start + overlap_length - 1, next_align.strand)
                    resolved_read_aligns.append(new_align)

                    # # update next align's cords in resolved_read_aligns
                    next_align.ref_start += overlap_length - 1
                    next_align.read_start += overlap_length - 1

                    # print("partitally overlapped by major, change this to major align")
                    # # update major align
                    major_align = next_align

                    # print("ERROR: non overlap type")
                    # exit()

        # return resolved_read_aligns
        return sorted(resolved_read_aligns, key=lambda x: (x.read_start, x.read_end))


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
            ops = ops[1: ]
            lengths = lengths[1: ]

        if ops[-1] == "S" or ops[-1] == "H":
            ops = ops[: -1]
            lengths = lengths[: -1]

    return ops, lengths
