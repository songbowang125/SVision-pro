import numpy as np

from src.cigar_op import collect_cigars_from_bam, collect_cigar_from_ref


class CigarMap:
    """
    cigar map for a ref position.
    self.map includes the detail cigars that contribute to this map
    self.cnt_map includes the number (support) of detail cigars in this map
    self.freq_map includes freq of self.cnt_map
    """
    def __init__(self):
        self.map = {"A": [], "T": [], "C": [], "G": [], "I": [], "D": [], "V": []}
        self.cnt_map = {"A": [], "T": [], "C": [], "G": [], "I": [], "D": [], "V": []}
        self.freq_map = {"A": 0, "T": 0, "C": 0, "G": 0, "I": 0, "D": 0, "V": 0}

    def clear_op(self, op):
        self.map[op] = None
        self.cnt_map[op] = 0
        self.freq_map[op] = 0

    def get_most_frequent_op(self):
        """
        given the freq max, return the most frequent op
        :param target_freq_map:
        :return:
        """
        most_freq_val = -1
        most_freq_op = ""

        for op in self.freq_map:
            if self.freq_map[op] > most_freq_val:
                most_freq_op = op
                most_freq_val = self.freq_map[op]

        return most_freq_op

    def sort_by_freq(self):
        """
        sort ops by freq
        :return:
        """

        non_zero_ops = []

        remaining_freq = 100

        # # extract non-zero ops
        for op in self.freq_map:
            op_freq = self.freq_map[op]

            if op_freq > 0:
                remaining_freq -= op_freq
                non_zero_ops.append([op, op_freq])

        non_zero_ops = sorted(non_zero_ops, key=lambda x: x[1], reverse=True)

        # # not homo, then we need add ref
        if remaining_freq > 0:
            non_zero_ops.append(["REF", remaining_freq])

        return non_zero_ops

    def to_string(self):
        return self.freq_map


def generate_cigar_map_from_ref(non_match_cigars, map_ref_chrom, map_ref_start, map_ref_end, options):
    """
    mapping cigar to map (matrix) for ref
    :param non_match_cigars:
    :param map_ref_chrom:
    :param map_ref_start:
    :param map_ref_end:
    :param options:
    :return:
    """
    cigar_maps = [CigarMap() for i in range(map_ref_start, map_ref_end + 1)]

    # # STEP: generate cigar map cnt
    for cigar in non_match_cigars:

        if map_ref_start <= cigar.ref_start <= cigar.ref_end <= map_ref_end:
            index = cigar.ref_start - map_ref_start
            cigar_maps[index].freq_map[cigar.op] = 100

    return cigar_maps


def generate_cigar_map_from_bam(hyper_cigars, map_ref_chrom, map_ref_start, map_ref_end, map_type, options):
    """
    mapping cigar to map (matrix)
    :return:
    """
    # # STEP: generate coverage map by detail cigar's coverage
    coverage_maps = np.ones(map_ref_end - map_ref_start + 1)
    for hyper_cigar in hyper_cigars:
        for detail_cigar in hyper_cigar.detail_cigars:
            if map_type == "base" and (detail_cigar.ref_end > map_ref_end or detail_cigar.ref_start < map_ref_start):
                continue
            
            coverage_maps[detail_cigar.ref_start - map_ref_start: detail_cigar.ref_end - map_ref_start + 1] = detail_cigar.coverage

            if map_type == "base" and detail_cigar.op == "MappedI":
                if detail_cigar.secondary_ref_end > map_ref_end or detail_cigar.secondary_ref_start < map_ref_start:
                    continue
                coverage_maps[detail_cigar.secondary_ref_start - map_ref_start: detail_cigar.secondary_ref_end - map_ref_start + 1] = detail_cigar.coverage
    
    # # generate cigar map
    cigar_maps = [CigarMap() for i in range(map_ref_start, map_ref_end + 1)]

    # # STEP: generate cigar map cnt
    for hyper_cigar in hyper_cigars:
        for detail_cigar in hyper_cigar.detail_cigars:
            if map_ref_start <= detail_cigar.ref_start <= detail_cigar.ref_end <= map_ref_end:
                cigar_op = "I" if detail_cigar.op in ["MappedI", "UnmappedI"] else detail_cigar.op

                index_start = detail_cigar.ref_start - map_ref_start
                index_end = detail_cigar.ref_end - map_ref_start

                # # increase supp num for each postion
                for i in range(index_start, index_end + 1):
                    # # update cnt map
                    cigar_maps[i].cnt_map[cigar_op].extend(hyper_cigar.supp_reads)

                    # # update map
                    if map_type == "base" and detail_cigar not in cigar_maps[i].map[cigar_op]:
                        cigar_maps[i].map[cigar_op].append(detail_cigar)
    
                # # for mapped I, we need mark the source positions
                if map_type == "base" and detail_cigar.op == "MappedI":
                    # # the secondary region is out of map's region, then continue

                    if detail_cigar.secondary_ref_end < map_ref_start or detail_cigar.secondary_ref_start > map_ref_end:
                        continue
                    
                    secondary_index_start = max(0, detail_cigar.secondary_ref_start - map_ref_start)
                    secondary_index_end = min(map_ref_end - map_ref_start, detail_cigar.secondary_ref_end - map_ref_start)

                    for i in range(secondary_index_start, secondary_index_end + 1):
                        cigar_maps[i].map[cigar_op].append(detail_cigar)
                        cigar_maps[i].cnt_map[cigar_op].extend(hyper_cigar.supp_reads)

    # # STEP: generate cigar map cnt
    for pos in range(len(cigar_maps)):
        for op in ["A", "T", "C", "G", "I", "D", "V"]:
            cigar_maps[pos].cnt_map[op] = len(set(cigar_maps[pos].cnt_map[op]))

    # # STEP: generate cigar map freq
    for pos in range(len(cigar_maps)):
        tot_read_num = coverage_maps[pos]
        # if map_type == "base":
        #     print(pos, tot_read_num, cigar_maps[pos].cnt_map["I"])
        if tot_read_num == 0:
            tot_read_num = sum([cigar_maps[pos].cnt_map[op] for op in cigar_maps[pos].cnt_map.keys()])
        if tot_read_num == 0:
            continue

        # # calculate freq
        for op in ["A", "T", "C", "G", "I", "D", "V"]:
            cigar_maps[pos].freq_map[op] = round(cigar_maps[pos].cnt_map[op] / tot_read_num * 100)

            # # NOTE: this aims to solve the small dels due to the imprecise sequencing, since they will make the coverage smaller than real
            if cigar_maps[pos].freq_map[op] > 100:
                cigar_maps[pos].freq_map[op] = 100

    return cigar_maps


def generate_cigar_map_for_partition(partition, options):
    """
    generate both base and target cigar map for this parittion
    :return:
    """

    # # STEP: set map's ref start and end to partition's
    map_ref_chrom, map_ref_start, map_ref_end = partition.ref_chrom, partition.left_extension_start, partition.right_extension_end

    # # STEP: collect from ref
    ref_hyper_cigars = collect_cigar_from_ref(options.genome_path, partition.ref_chrom, partition.left_extension_start, partition.right_extension_end, options)
    ref_cigar_maps = generate_cigar_map_from_ref(ref_hyper_cigars, map_ref_chrom, map_ref_start, map_ref_end, options)

    # # STEP: generate targe map
    target_cigar_maps = generate_cigar_map_from_bam(partition.included_hyper_cigars, map_ref_chrom, map_ref_start, map_ref_end, "target", options)

    # # STEP: generate base map
    base_cigar_maps = {}

    if len(options.base_path) != 0:
        for base_file_index in partition.base_partitions.keys():

            if len(partition.base_partitions[base_file_index]) == 0:
                base_cigar_maps[base_file_index] = ref_cigar_maps

            else:
                # # STEP: collect merged hyper cigar in each partition
                merged_base_hyper_cigars = []
                for base_partition in partition.base_partitions[base_file_index]:
                    merged_base_hyper_cigars.extend(base_partition.included_hyper_cigars)

                # # STEP: generate cigar map for base
                base_cigar_maps[base_file_index] = generate_cigar_map_from_bam(merged_base_hyper_cigars, map_ref_chrom, map_ref_start, map_ref_end, "base", options)

    else:
        base_cigar_maps[0] = ref_cigar_maps

    # # STEP: set
    partition.set_cigar_maps(target_cigar_maps, ref_cigar_maps, base_cigar_maps)


def generate_gt_cigar_map_for_partition(partition, cigars_gt, options):
    """
    pass
    :param var_gt:
    :return:
    """

    # # STEP: set map's ref start and end to partition's
    map_ref_chrom, map_ref_start, map_ref_end = partition.ref_chrom, partition.left_extension_start, partition.right_extension_end

    base_cigar_maps = partition.base_cigar_maps

    gt_target_cigar_maps = [CigarMap() for i in range(map_ref_start, map_ref_end + 1)]
    gt_base_cigar_maps = [CigarMap() for i in range(map_ref_start, map_ref_end + 1)]

    # # STEP: generate cigar map
    for cur_gt in cigars_gt:
        gt_cigar_op = cur_gt[0]
        gt_ref_start = cur_gt[1] - map_ref_start
        gt_ref_end = cur_gt[2] - map_ref_start
        gt_var_type = cur_gt[3]

        for index in range(gt_ref_start, gt_ref_end + 1):
            # # STEP: set target map
            gt_target_cigar_maps[index].cnt_map[gt_cigar_op] = 100
            gt_target_cigar_maps[index].freq_map[gt_cigar_op] = 100

            # # STEP: set base map accroding to germline or somatic
            if gt_var_type == "somatic":
                # # somatic means gt base map is same as base map
                gt_base_cigar_maps[index] = base_cigar_maps[index]

            elif gt_var_type == "germline":
                # # germline means gt base map is same as gt target map
                gt_base_cigar_maps[index] = gt_target_cigar_maps[index]

            else:
                pass

    # # STEP: set
    partition.set_gt_cigar_maps(gt_target_cigar_maps, gt_base_cigar_maps)

