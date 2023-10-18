import collections
import pysam
from src.cigar_op import HyperCigar, DetailCigar
import numpy as np
from src.cluster_op import hierarchical_cluster, iterative_cluster, dbscan_cluster
from src.align_op import fetch_ref_seq


class Partition:
    def __init__(self, included_hyper_cigars):

        self.update_included_hyper_cigars(included_hyper_cigars)

        self.bonus_supp_reads = []      # bonus reads include: uncovered and shorter supp reads
        self.filter_marker = "PASS"
        self.right_extension_end = None
        self.left_extension_start = None
        self.base_cigar_maps = None
        self.ref_cigar_maps = None
        self.target_cigar_maps = None
        self.gt_base_cigar_maps = None
        self.gt_target_cigar_maps = None
        self.target_variant = []
        self.target_variant_type = "NA"
        self.base_variants = {}
        self.base_variants_types = "NA"
        self.base_partitions = {}

        self.color_plot_resize_ratio = 1
        self.color_plot_name = ""
        self.genotype = "./."

    def set_filter_marker(self, marker):
        self.filter_marker = marker

    def merge_included_hyper_cigars(self):
        """
        use voting strategy to merge separate hyper cigars to final one
        """

        included_hyper_cigar_ops = []
        for hyper_cigar in self.included_hyper_cigars:
            included_hyper_cigar_ops.extend([hyper_cigar.op for i in range(len(hyper_cigar.supp_reads))])

        # # a mixed partition, such as I, I+I+I were merged together, then we only kept the most common one
        if len(set(included_hyper_cigar_ops)) != 1:

            most_common_cigar_ops = collections.Counter(included_hyper_cigar_ops).most_common()[0][0]

            origin_supp_reads = self.supp_reads
            new_included_hyper_cigars = []

            # # only keep the most common, but keep all supp reads
            for hyper_cigar in self.included_hyper_cigars:
                if hyper_cigar.op == most_common_cigar_ops:
                    new_included_hyper_cigars.append(hyper_cigar)

            self.update_included_hyper_cigars(new_included_hyper_cigars)

            self.supp_reads = origin_supp_reads

        # # STEP: begin merge
        best_detail_cigars = []

        # # STEP: deal with each op (detail cigars)
        nested_op_split = self.included_hyper_op.split("+")
        for i in range(len(nested_op_split)):

            all_detail_ops = []
            all_ref_chroms = []
            all_ref_starts = []
            all_ref_ends = []
            all_hybrid_lengths = []

            all_second_ref_chroms = []
            all_second_ref_starts = []
            all_second_ref_ends = []
            all_second_ref_strands = []

            all_bnd_ref_chroms = []
            all_bnd_ref_starts = []
            all_bnd_ref_ends = []
            all_bnd_ref_strands = []

            all_alt_bases = []
            # # collect all corresponding infos
            for hyper_cigar in self.included_hyper_cigars:
                # print(hyper_cigar.to_string(), i)
                detail_cigar = hyper_cigar.detail_cigars[i]

                # # each supp read means an event
                for j in range(len(hyper_cigar.supp_reads)):
                    all_alt_bases.append(detail_cigar.alt_bases)

                    all_detail_ops.append(detail_cigar.op)
                    all_ref_chroms.append(detail_cigar.ref_chrom)
                    all_ref_starts.append(detail_cigar.ref_start)
                    all_ref_ends.append(detail_cigar.ref_end)
                    all_hybrid_lengths.append(detail_cigar.hybrid_length)

                    all_second_ref_chroms.append(detail_cigar.secondary_ref_chrom)
                    all_second_ref_starts.append(detail_cigar.secondary_ref_start)
                    all_second_ref_ends.append(detail_cigar.secondary_ref_end)
                    all_second_ref_strands.append(detail_cigar.secondary_ref_strand)

                    all_bnd_ref_chroms.append(detail_cigar.bnd_ref_chrom)
                    all_bnd_ref_starts.append(detail_cigar.bnd_ref_start)
                    all_bnd_ref_ends.append(detail_cigar.bnd_ref_end)
                    all_bnd_ref_strands.append(detail_cigar.bnd_ref_strand)

            # # choose the most common one as the refined info
            best_detail_op = max(set(all_detail_ops), key=all_detail_ops.count)
            best_ref_chrom = max(set(all_ref_chroms), key=all_ref_chroms.count)

            all_ref_starts = sorted(all_ref_starts)
            ref_start_stats = collections.Counter(all_ref_starts).most_common()
            if ref_start_stats[0][1] >= 5:
                best_ref_start = ref_start_stats[0][0]
            else:
                best_ref_start = all_ref_starts[int(len(all_ref_starts) / 2)]

            all_hybrid_lengths = sorted(all_hybrid_lengths)

            hybrid_length_stats = collections.Counter(all_hybrid_lengths).most_common()
            if hybrid_length_stats[0][1] >= 5:
                best_hybrid_length = hybrid_length_stats[0][0]
            else:
                best_hybrid_length = all_hybrid_lengths[int(len(all_hybrid_lengths) / 2)]

            # print(self.ref_start, self.hybrid_length, len(self.supp_reads))
            # for i in range(len(all_alt_bases)):
            #     print(len(all_alt_bases[i]), all_hybrid_lengths[i])

            best_alt_bases = all_alt_bases[all_hybrid_lengths.index(best_hybrid_length)]

            # # newly add in v0.11.20, refine the length
            if best_detail_op not in ["MappedI", "UnmappedI", "B"]:
                best_ref_end = best_ref_start + best_hybrid_length
            else:
                best_ref_end = best_ref_start
            # # in case that: None are the most frequent cords for mappedI
            if best_detail_op == "MappedI":
                all_second_ref_chroms = [val for val in all_second_ref_chroms if val is not None]
                all_second_ref_starts = [val for val in all_second_ref_starts if val is not None]
                all_second_ref_ends = [val for val in all_second_ref_ends if val is not None]
                all_second_ref_strands = [val for val in all_second_ref_strands if val is not None]

            if best_detail_op == "B":
                all_bnd_ref_chroms = [val for val in all_bnd_ref_chroms if val is not None]
                all_bnd_ref_starts = [val for val in all_bnd_ref_starts if val is not None]
                all_bnd_ref_ends = [val for val in all_bnd_ref_ends if val is not None]
                all_bnd_ref_strands = [val for val in all_bnd_ref_strands if val is not None]

            best_second_ref_chrom = max(set(all_second_ref_chroms), key=all_second_ref_chroms.count)
            best_second_ref_start = max(set(all_second_ref_starts), key=all_second_ref_starts.count)
            best_second_ref_end = max(set(all_second_ref_ends), key=all_second_ref_ends.count)
            best_second_ref_strand = max(set(all_second_ref_strands), key=all_second_ref_strands.count)

            best_bnd_ref_chrom = max(set(all_bnd_ref_chroms), key=all_bnd_ref_chroms.count)
            best_bnd_ref_start = max(set(all_bnd_ref_starts), key=all_bnd_ref_starts.count)
            best_bnd_ref_end = max(set(all_bnd_ref_ends), key=all_bnd_ref_ends.count)
            best_bnd_ref_strand = max(set(all_bnd_ref_strands), key=all_bnd_ref_strands.count)

            # # generate best detail cigar
            tmp_detail_cigar = DetailCigar(best_detail_op, best_hybrid_length, best_ref_chrom, best_ref_start, best_ref_end, "+")
            tmp_detail_cigar.set_alt_bases(best_alt_bases)
            tmp_detail_cigar.set_secondary_mapping(best_second_ref_chrom, best_second_ref_start, best_second_ref_end, best_second_ref_strand)
            tmp_detail_cigar.set_bnd_mapping(best_bnd_ref_chrom, best_bnd_ref_start, best_bnd_ref_end, best_bnd_ref_strand)

            best_detail_cigars.append(tmp_detail_cigar)

        # # update hyper cigar
        self.update_included_hyper_cigars([HyperCigar(best_detail_cigars, self.supp_reads)])

    def add_included_hyper_cigars(self, new_hyper_cigars):
        hyper_cigars = self.included_hyper_cigars + new_hyper_cigars

        self.update_included_hyper_cigars(hyper_cigars)

    def update_included_hyper_cigars(self, included_hyper_cigars):

        hyper_cigar_num = len(included_hyper_cigars)
        medine_index = int(hyper_cigar_num / 2)

        self.included_hyper_cigars = included_hyper_cigars
        self.included_hyper_op = included_hyper_cigars[medine_index].op

        self.ref_chrom = included_hyper_cigars[0].ref_chrom
        self.ref_start = sorted([cigar.ref_start for cigar in included_hyper_cigars])[medine_index]
        self.ref_end = sorted([cigar.ref_end for cigar in included_hyper_cigars])[medine_index]

        self.ref_start_with_secondary = sorted([cigar.ref_start_with_secondary for cigar in included_hyper_cigars])[medine_index]
        self.ref_end_with_secondary = sorted([cigar.ref_end_with_secondary for cigar in included_hyper_cigars])[medine_index]

        self.hybrid_length = sorted([cigar.hybrid_length for cigar in included_hyper_cigars])[medine_index]
        self.read_length = self.calculate_read_length()
        self.ref_length = self.ref_end - self.ref_start

        self.supp_reads = []
        for cigar in self.included_hyper_cigars:
            self.supp_reads.extend(cigar.supp_reads)
        self.supp_reads = list(set(self.supp_reads))

    def set_base_partitions(self, base_path_index, base_partitions):
        
        self.base_partitions[base_path_index] = base_partitions
        
        # # deal with near Is but are mapped in different position (very close ) due to mapping artifact
        self_p_op_unique = set(self.included_hyper_cigars[0].op.split("+"))
        if len(self_p_op_unique) == 1 and list(self_p_op_unique)[0] == "I":
            for base_p in self.base_partitions[base_path_index]:
                # # near in 10 bp
                base_p_op_unique = set(base_p.included_hyper_cigars[0].op.split("+"))

                if len(base_p_op_unique) == 1 and list(base_p_op_unique)[0] == "I" and abs(self.included_hyper_cigars[0].ref_start - base_p.included_hyper_cigars[0].ref_start) < 500:
                    if self.included_hyper_cigars[0].hybrid_length <= base_p.included_hyper_cigars[0].hybrid_length and (self.included_hyper_cigars[0].hybrid_length / base_p.included_hyper_cigars[0].hybrid_length) >= 0.8:
                        # only update inserted position
                        base_p.included_hyper_cigars[0].ref_start = self.included_hyper_cigars[0].ref_start
                        base_p.included_hyper_cigars[0].ref_end = self.included_hyper_cigars[0].ref_end
                        base_p.included_hyper_cigars[0].detail_cigars[0].ref_start = self.included_hyper_cigars[0].detail_cigars[0].ref_start
                        base_p.included_hyper_cigars[0].detail_cigars[0].ref_end = self.included_hyper_cigars[0].detail_cigars[0].ref_end
                    if self.included_hyper_cigars[0].hybrid_length > base_p.included_hyper_cigars[0].hybrid_length and (base_p.included_hyper_cigars[0].hybrid_length / self.included_hyper_cigars[0].hybrid_length) >= 0.8:
                        # only update inserted position
                        base_p.included_hyper_cigars[0].ref_start = self.included_hyper_cigars[0].ref_start
                        base_p.included_hyper_cigars[0].ref_end = self.included_hyper_cigars[0].ref_end
                        base_p.included_hyper_cigars[0].detail_cigars[0].ref_start = self.included_hyper_cigars[0].detail_cigars[0].ref_start
                        base_p.included_hyper_cigars[0].detail_cigars[0].ref_end = self.included_hyper_cigars[0].detail_cigars[0].ref_end

    def calculate_read_length(self):
        read_length = 0
        for cigar in self.included_hyper_cigars:
            # # only D dose not need to increase read length
            if cigar.op != "D":
                read_length += cigar.hybrid_length

        return read_length

    def update_alt_bases(self, ref_file):
        self.included_hyper_cigars[0].alt_bases = ""

        for detail_cigar in self.included_hyper_cigars[0].detail_cigars:
            if detail_cigar.alt_bases == "":
                detail_cigar.alt_bases = fetch_ref_seq(detail_cigar.ref_chrom, detail_cigar.ref_start, detail_cigar.ref_end, ref_file)

            self.included_hyper_cigars[0].alt_bases += detail_cigar.alt_bases

    def set_cigar_maps(self, target_cigar_maps, ref_cigar_maps, base_cigar_maps):
        self.target_cigar_maps = target_cigar_maps
        self.ref_cigar_maps = ref_cigar_maps
        self.base_cigar_maps = base_cigar_maps

    def set_gt_cigar_maps(self, gt_target_cigar_maps, gt_base_cigar_maps):
        self.gt_target_cigar_maps = gt_target_cigar_maps
        self.gt_base_cigar_maps = gt_base_cigar_maps

    def set_target_variant(self, target_variant):
        self.target_variant = target_variant
        self.target_variant_type = "+".join([var.type for var in self.target_variant])

    def set_base_variants(self, base_index, base_variants):
        self.base_variants[base_index] = base_variants
        self.base_variants_types = {}

        for base_file_index in self.base_variants.keys():
            self.base_variants_types[base_file_index] = []
            for base_variant in self.base_variants[base_file_index]:
                self.base_variants_types[base_file_index].append("+".join([var.type for var in base_variant]))

    def set_extension(self, left_extension_start, right_extension_end):
        self.left_extension_start = left_extension_start
        self.right_extension_end = right_extension_end

    def set_color_plot_name(self, name):
        self.color_plot_name = name

    def set_color_plot_resize_ratio(self, resize_ratio):
        self.color_plot_resize_ratio = resize_ratio

    def get_detail_cigars(self):
        detail_cigars = []
        for hyper_cigar in self.included_hyper_cigars:
            detail_cigars.extend(hyper_cigar.detail_cigars)

        return detail_cigars

    def is_cover_new_hyper_cigar(self, new_hyper_cigar, min_sv_size):
        # # different op
        if self.included_hyper_op != new_hyper_cigar.op:
            return False
        elif self.ref_chrom != new_hyper_cigar.ref_chrom:
            return False
        elif abs(self.hybrid_length - new_hyper_cigar.hybrid_length) > min_sv_size:
            return False
        # # same op, then we need consider I and others individually
        else:
            # # process I
            if new_hyper_cigar.op == "I":
                # # for I, we need make sure the detail cigar ops (MappedI or unmappedI) are same
                if new_hyper_cigar.detail_cigars[0].op == self.included_hyper_cigars[0].detail_cigars[0].op and new_hyper_cigar.ref_start - (self.ref_end + self.hybrid_length) <= 1:
                    return True
            elif new_hyper_cigar.op == "B":
                if (new_hyper_cigar.ref_start - self.ref_start <= 1000 or new_hyper_cigar.detail_cigars[0].bnd_ref_start - self.included_hyper_cigars[0].detail_cigars[0].bnd_ref_start <= 1000) and new_hyper_cigar.detail_cigars[0].bnd_ref_chrom == self.included_hyper_cigars[0].detail_cigars[0].bnd_ref_chrom and new_hyper_cigar.op == self.included_hyper_op:
                    return True
            else:
                if new_hyper_cigar.ref_start - self.ref_end <= 10 and new_hyper_cigar.op == self.included_hyper_op:
                    return True

        return False

    def to_string(self):
        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.ref_chrom, self.ref_start, self.ref_end, len(self.included_hyper_cigars), [cigar.to_string() for cigar in self.included_hyper_cigars], self.cigar_map_freq)


def perform_best_partition_cluster(fine_partitions, method="hierarchical"):
    """
    fine partition, use perform hierarchical clustering
    """


    candidate_partitions = []

    for fine_partition in fine_partitions:
        included_hyper_cigars = fine_partition.included_hyper_cigars

        if len(included_hyper_cigars) == 1:
            candidate_partitions.append(Partition(included_hyper_cigars))

        else:
            # if coarse_partition == "B":
            #     cluster_data = np.array([[hyper_cigar.ref_start, hyper_cigar.ref_start + 1000, 1000] for hyper_cigar in included_hyper_cigars])
            # else:
            cluster_data = np.array([[hyper_cigar.ref_start, hyper_cigar.ref_start + hyper_cigar.hybrid_length, 1000] for hyper_cigar in included_hyper_cigars])

            if method == "hierarchical":
                cluster_indices = hierarchical_cluster(cluster_data, 0.4, )
                # # STEP: parse res from hierarchical cluster
                new_clusters = [[] for i in range(max(cluster_indices))]
                for cigar_index, cluster_index in enumerate(cluster_indices):
                    new_clusters[cluster_index - 1].append(included_hyper_cigars[cigar_index])

                for hyper_cigars_in_this_cluster in new_clusters:
                    candidate_partitions.append(Partition(hyper_cigars_in_this_cluster))

            elif method == "dbscan":
                pass
                # cluster_indices = dbscan_cluster(cluster_data)
                #
                # # # STEP: parse res from dbscan cluster
                # new_clusters = [[] for i in range(max(cluster_indices) + 1)]
                #
                # for cigar_index, cluster_index in enumerate(cluster_indices):
                #     new_clusters[cluster_index].append(included_hyper_cigars[cigar_index])
                #
                # for hyper_cigars_in_this_cluster in new_clusters:
                #     candidate_partitions.append(Partition(hyper_cigars_in_this_cluster))

            else:
                pass

    return candidate_partitions


def perform_coarse_partition_iter(hyper_cigars, options):
    """
    use iterative method to cluster, the coarse partition does not care op types
    """

    coarse_partitions = []

    uncovered_hyper_cigars = []
    shorter_hyper_cigars = []

    if len(hyper_cigars) == 0:
        return coarse_partitions

    #hyper_cigars = sorted(hyper_cigars, key=lambda x: x.hybrid_length)
    hyper_cigars = sorted(hyper_cigars, key=lambda x: len(x.supp_reads), reverse=True)

    # # STEP: traverse all non match cigars and find in-continuous cigar
    for cur_cigar in hyper_cigars:

        # # leave uncovered cigar alone for later process
        if not options.rescue_large and cur_cigar.uncovered_flag is True:
            uncovered_hyper_cigars.append(cur_cigar)

        elif options.rescue_large and cur_cigar.uncovered_flag is True and len(cur_cigar.supp_reads) < 3:
            uncovered_hyper_cigars.append(cur_cigar)

        elif cur_cigar.shorter_flag is True:
            shorter_hyper_cigars.append(cur_cigar)
        else:
            candidate_merge = []
            for i in range(len(coarse_partitions)):
                partition = coarse_partitions[i]
                merge_flag, merge_sim = iterative_cluster(partition, cur_cigar, options, ignore_op=True)
                if merge_flag is True:
                    candidate_merge.append([i, merge_sim])

            # # if covered, add this cigar to partition
            if len(candidate_merge) != 0:
                best_merge_partition_index = sorted(candidate_merge, key=lambda x: x[1], reverse=True)[0][0]

                coarse_partitions[best_merge_partition_index].add_included_hyper_cigars([cur_cigar])

            # # if not covered, generate a new partition
            else:
                coarse_partitions.append(Partition([cur_cigar]))

    coarse_partitions = sorted(coarse_partitions, key=lambda x: len(x.supp_reads), reverse=True)

    # # STEP: for each shorter cigar, merge it directly into partition
    for shorter_hyper_cigar in shorter_hyper_cigars:
        for coarse_partition in coarse_partitions:
            if len(coarse_partition.supp_reads) < options.min_supp:
                continue
            if shorter_hyper_cigar.op == "I" and coarse_partition.included_hyper_op == shorter_hyper_cigar.op and abs(coarse_partition.ref_start - shorter_hyper_cigar.ref_start) < 20 and abs(coarse_partition.hybrid_length - shorter_hyper_cigar.hybrid_length) < 20:
                coarse_partition.supp_reads.extend(shorter_hyper_cigar.supp_reads)
                coarse_partition.bonus_supp_reads.extend(shorter_hyper_cigar.supp_reads)

                break
            elif shorter_hyper_cigar.op == "D" and coarse_partition.included_hyper_op == shorter_hyper_cigar.op and abs(coarse_partition.ref_start - shorter_hyper_cigar.ref_start) < 50 and abs(coarse_partition.hybrid_length - shorter_hyper_cigar.hybrid_length) < 20:
                coarse_partition.supp_reads.extend(shorter_hyper_cigar.supp_reads)
                coarse_partition.bonus_supp_reads.extend(shorter_hyper_cigar.supp_reads)

                break

    if options.detect_mode == "germline" or options.detect_mode == 'genotype' or options.force_cluster:
        # # STEP: for each uncovered cigar, merge it directly into partition

        for uncovered_hyper_cigar in uncovered_hyper_cigars:
            for coarse_partition in coarse_partitions:
                if "I" in coarse_partition.included_hyper_op:
                    merge_flag, merge_sim = iterative_cluster(coarse_partition, uncovered_hyper_cigar, options, ignore_op=True, ignore_length=True)

                    if merge_flag is True:
                        # print("merge into", coarse_partition.included_hyper_cigars[0].hybrid_length)
                        coarse_partition.supp_reads.extend(uncovered_hyper_cigar.supp_reads)
                        coarse_partition.bonus_supp_reads.extend(uncovered_hyper_cigar.supp_reads)

                        break

    return coarse_partitions


def perform_fine_partition_iter(coarse_partitions, options):
    """
    split coarse partitions and get the fine partition
    """

    fine_partitions = []

    for coarse_partition in coarse_partitions:

        included_hyper_cigars = coarse_partition.included_hyper_cigars
        included_hyper_cigar_ops = [hyper_cigar.op for hyper_cigar in included_hyper_cigars]
        bonus_supp_reads = coarse_partition.bonus_supp_reads

        # # STEP: simplify cigar op, or remove redundant sub-ops, e.g. I+I+I+V --> I+V
        included_hyper_cigar_ops_simple = []
        for hyper_cigar_op in included_hyper_cigar_ops:
            hyper_cigar_op = hyper_cigar_op.split("+")

            hyper_cigar_op_simple = list(set(hyper_cigar_op))
            hyper_cigar_op_simple.sort(key=hyper_cigar_op.index)

            included_hyper_cigar_ops_simple.append("+".join(hyper_cigar_op_simple))

        # # STEP: if there is only one kind of cigar op after simplification, there this is a fine partition
        # # otherwise, we need to split this coarse partition
        if len(set(included_hyper_cigar_ops_simple)) == 1:
            fine_partitions.append(coarse_partition)
        else:
            splited_partitions = []

            # # more than one op, we need split the origin partition
            for cigar_op_simple in list(set(included_hyper_cigar_ops_simple)):
                included_hyper_cigars_split = []
                for i in range(len(included_hyper_cigar_ops_simple)):
                    if included_hyper_cigar_ops_simple[i] == cigar_op_simple:
                        included_hyper_cigars_split.append(included_hyper_cigars[i])

                splited_partitions.append(Partition(included_hyper_cigars_split))

            splited_partitions = sorted(splited_partitions, key=lambda x: len(x.supp_reads), reverse=True)
            # # add the uncovered supp reads to the most supported one
            if len(splited_partitions[0].supp_reads) >= options.min_supp:
                splited_partitions[0].supp_reads.extend(bonus_supp_reads)

            fine_partitions.extend(splited_partitions)

    return fine_partitions


def generate_partitions_in_interval(target_hyper_cigars, ref_file, interval_chrom, interval_start, interval_end, mode, options):
    """
    generate partitions from hyper cigars, there are two steps:
    1. according to their ops, split them (coarse partition)
    2. in each op type, perform hierarchical clustering (fine partition)
    :return:
    """

    # # STEP: perform partition
    coarse_partitions = perform_coarse_partition_iter(target_hyper_cigars, options)

    fine_partitions = perform_fine_partition_iter(coarse_partitions, options)

    # best_partitions = perform_best_partition_cluster(fine_partitions)

    # # STEP: filter    
    filtered_partitions = []
    for partition in fine_partitions:
        partition_info = "{}_{}_{}_{}, {}_{}, {}".format(partition.ref_chrom, partition.ref_start, partition.ref_end, partition.hybrid_length, partition.left_extension_start, partition.right_extension_end, partition.hybrid_length)

        # # STEP: merge include hyper cigar list to final one
        partition.merge_included_hyper_cigars()

        if partition.ref_start > partition.ref_end:
            continue
        if partition.ref_end - partition.ref_start > options.max_sv_size:
            continue

        # # STEP: filter
        if partition.included_hyper_op == "B":
            if options.min_supp <= len(partition.supp_reads):
                filtered_partitions.append(partition)
        else:
            if mode == "target" and options.min_sv_size <= partition.hybrid_length <= options.max_sv_size and options.min_supp <= len(partition.supp_reads):
                partition.update_alt_bases(ref_file)
                if "N" not in partition.included_hyper_cigars[0].alt_bases and "N" not in ref_file.fetch(partition.ref_chrom, max(partition.ref_start - options.dist_diff_length, 0), partition.ref_start + options.dist_diff_length):
                    filtered_partitions.append(partition)

            if mode == "base" and options.min_sv_size <= partition.hybrid_length <= options.max_sv_size and 1 <= len(partition.supp_reads):
                partition.update_alt_bases(ref_file)
                if "N" not in partition.included_hyper_cigars[0].alt_bases and "N" not in ref_file.fetch(partition.ref_chrom, max(partition.ref_start - options.dist_diff_length, 0), partition.ref_start + options.dist_diff_length):
                    filtered_partitions.append(partition)

    # # STEP: mark imprecise partition
    mark_poly_partitions(filtered_partitions, interval_chrom, interval_start, interval_end)

    # # STEP: extend partition
    filtered_partitions = extend_partitions(filtered_partitions)

    return filtered_partitions


def extend_partitions(candidate_partitions):
    """
    extend signature, mainly searching matching extensions at both the left and right sides
    :param candidate_partitions:
    :param options:
    :return:
    """

    candidate_partitions = sorted(candidate_partitions, key=lambda x: x.ref_start)
    for i in range(len(candidate_partitions)):

        cur_partition = candidate_partitions[i]

        # # set max extension same as patition length, for complex/nested event, we extend less for clear representation
        # if "+" in cur_partition.included_hyper_op:
        max_extension = int(0.7 * (cur_partition.hybrid_length + 1))
        # else:
        #     max_extension = int(1.4 * (cur_partition.hybrid_length + 1))

        left_extension = min(cur_partition.ref_start - max_extension, cur_partition.ref_start_with_secondary)
        right_extension = max(cur_partition.ref_end + max_extension, cur_partition.ref_end_with_secondary)
        # left_extension = cur_partition.ref_start - max_extension
        # right_extension = cur_partition.ref_end + max_extension

        # # fix when the left or right extension length are so long
        left_extension_length = abs(cur_partition.ref_start - left_extension)
        right_extension_length = abs(right_extension - cur_partition.ref_end)

        # left extension length is too long, extend right more
        extension_ratio = left_extension_length / right_extension_length
        if extension_ratio > 5:
            right_extension = right_extension + int(round(extension_ratio / 5) * right_extension_length)

        # right extension length is too long, extend left more
        extension_ratio = right_extension_length / left_extension_length
        if extension_ratio > 5:
            left_extension = left_extension - int(round(extension_ratio / 5) * left_extension_length)

        cur_partition.set_extension(left_extension, right_extension)

    return candidate_partitions


def mark_poly_partitions(partition_list, interval_chrom, interval_start, interval_end):
    """
    mark imprecise partitions. There are two imprecise conditions:
    1. BND
    2. more than two partition located in a small bin (e.g. 1000)
    """
    
    bin_size = 200

    bin_num = int((interval_end - interval_start) / bin_size) + 1
    
    # # STEP: save partition_index into bins
    bin_records = [[] for i in range(bin_num + 1)]

    for partition_index in range(len(partition_list)):
        if partition_list[partition_index].ref_chrom != interval_chrom:
            continue

        if not (interval_start <= partition_list[partition_index].ref_start <= partition_list[partition_index].ref_end <= interval_end):
            continue

        if partition_list[partition_index].included_hyper_op == "B":
            partition_list[partition_index].set_filter_marker("POLY")
        else:
            belonged_bin_index = int((partition_list[partition_index].ref_start - interval_start) / bin_size)

            bin_records[belonged_bin_index].append(partition_index)

    # # STEP: traverse bins to find POLY
    for bin in bin_records:
        if len(bin) > 1:
            for partition_index in bin:
                partition_list[partition_index].set_filter_marker("POLY")
