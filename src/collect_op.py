from src.cigar_op import collect_cigars_from_bam
from src.cigar_map import generate_cigar_map_for_partition
from src.hash_realign.run_hash_realign import hash_realign_unmapped
from src.partition_op import generate_partitions_in_interval
from src.variant_op import detect_target_pseudo_variant_for_partition, detect_base_pseudo_variant_for_partition
from src.color_plot import generate_color_plot_for_partition
import datetime
import os
from src.output_vcf import output_vcf_header_in_interval, output_origin_record_to_vcf
import traceback
import sys
import logging
import pysam


def collect_and_detect_in_interval(interval_chrom, interval_start, interval_end, options):
    """
    process each single partial interval
    """

    interval_start_time = datetime.datetime.now()

    ref_file = pysam.FastaFile(options.genome_path)

    # # STEP: generate partial file
    interval_out_path = os.path.join(options.sample_out_path, "{}_{}_{}".format(interval_chrom, interval_start, interval_end))
    if not os.path.exists(interval_out_path):
        os.mkdir(interval_out_path)
    interval_img_path = os.path.join(interval_out_path, "Images")
    if not os.path.exists(interval_img_path):
        os.mkdir(interval_img_path)

    interval_list_out = open(os.path.join(interval_out_path, "test.txt"), "w")
    interval_vcf_out = open(os.path.join(interval_out_path, "vcf.txt"), "w")

    output_vcf_header_in_interval(options.genome_path, options.target_path, options.base_path, interval_vcf_out)

    errored_interval = []
    errored_partition = []

    try:
        # # STEP: collect non match cigars from target file (tumor BAM)
        target_hyper_cigars, target_interval_coverage = collect_cigars_from_bam(options.target_path, interval_chrom, interval_start, interval_end, "target", options)

        # # STEP: run coarse grained segmentation
        partitions_list = generate_partitions_in_interval(target_hyper_cigars, ref_file, interval_chrom, interval_start, interval_end, "target", options)

    except:
        target_interval_coverage = []
        error_type, error_value, error_trace = sys.exc_info()
        error_log = "Error log: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
        errored_interval = ["{}_{}_{}".format(interval_chrom, interval_start, interval_end), error_log]

        return errored_interval, errored_partition

    # # STEP: for each coarse partition, do fine segmentation
    for partition in partitions_list:
        try:

            # if options.skip_inheritype:
            #     if "I+D" in partition.included_hyper_op:
            #         hash_realign_unmapped(partition, ref_file)

            # # STEP: detect variants
            detect_target_pseudo_variant_for_partition(partition, interval_start, interval_end, target_interval_coverage, options)

            if not options.skip_inheritype:

                for base_index in range(len(options.base_path)):
                    cur_base_path = options.base_path[base_index]

                    # # STEP: collect from base
                    if (partition.ref_end_with_secondary - partition.ref_start_with_secondary) < options.max_sv_size:
                        base_hyper_cigars, base_interval_coverage = collect_cigars_from_bam(cur_base_path, partition.ref_chrom, partition.left_extension_start, partition.right_extension_end, "base", options)
                    else:
                        base_hyper_cigars, base_interval_coverage = collect_cigars_from_bam(cur_base_path, partition.ref_chrom, partition.ref_start - 1000, partition.ref_end + 1000, "base", options)

                    # # STEP: generate partition for base hyper cigars
                    base_partitions = generate_partitions_in_interval(base_hyper_cigars, ref_file, partition.ref_chrom, partition.left_extension_start, partition.right_extension_end, "base", options)

                    target_hyper_cigar = partition.included_hyper_cigars[0]

                    base_partitions_similarities = []

                    for base_partition in base_partitions:
                        base_hyper_cigar = base_partition.included_hyper_cigars[0]

                        if set(base_hyper_cigar.op.split("+")) != set(target_hyper_cigar.op.split("+")):
                            op_similarity = 0
                        else:
                            op_similarity = 100

                        if op_similarity == 0:
                            base_partitions_similarities.append(0)
                            continue

                        if target_hyper_cigar.hybrid_length > base_hyper_cigar.hybrid_length:
                            size_similarity = round(100 * (base_hyper_cigar.hybrid_length / target_hyper_cigar.hybrid_length))
                        else:
                            size_similarity = round(100 * (target_hyper_cigar.hybrid_length / base_hyper_cigar.hybrid_length))

                        if size_similarity / 100 < options.size_sim_ratio - 0.1:
                            base_partitions_similarities.append(0)
                            continue

                        if target_hyper_cigar.hybrid_length < options.dist_diff_length:
                            target_range_length = options.dist_diff_length
                        else:
                            target_range_length = target_hyper_cigar.hybrid_length
                        target_range = set(range(target_hyper_cigar.ref_start, target_hyper_cigar.ref_start + target_range_length))

                        if base_hyper_cigar.hybrid_length < options.dist_diff_length:
                            base_range_length = options.dist_diff_length
                        else:
                            base_range_length = base_hyper_cigar.hybrid_length

                        base_range = range(base_hyper_cigar.ref_start, base_hyper_cigar.ref_start + base_range_length)
                        pos_similarity = int(100 * (len(target_range.intersection(base_range)) / target_range_length))

                        if pos_similarity == 0:
                            base_partitions_similarities.append(0)
                            continue

                        base_partitions_similarities.append(op_similarity + size_similarity + pos_similarity)

                    if len(base_partitions_similarities) == 0:
                        filtered_base_partitions = []
                    else:
                        max_similarity = max(base_partitions_similarities)

                        if max_similarity == 0:
                            filtered_base_partitions = []
                        else:
                            filtered_base_partitions = []

                            for i in range(len(base_partitions_similarities)):
                                if base_partitions_similarities[i] == max_similarity:
                                    filtered_base_partitions.append(base_partitions[i])

                    partition.set_base_partitions(base_index, filtered_base_partitions)

                    detect_base_pseudo_variant_for_partition(partition, partition.left_extension_start, partition.right_extension_end, base_index, base_interval_coverage, options)

                if (partition.ref_end_with_secondary - partition.ref_start_with_secondary) < options.max_sv_size:
                    # # STEP: generate base and target cigar maps
                    # # for bnd
                    if not options.skip_bnd and partition.included_hyper_cigars[0].op == "B":
                        pass
                    else:
                        generate_cigar_map_for_partition(partition, options)

                        # # STEP: generate color plots and write img name to file
                        color_plot_names = generate_color_plot_for_partition(partition, interval_img_path, options)
                        partition.set_color_plot_name(color_plot_names)
                        for plot_name in color_plot_names:
                            interval_list_out.write("{}\n".format(plot_name))

            # # STEP: output to vcf
            output_origin_record_to_vcf(partition, interval_vcf_out, options)

        except:
            error_type, error_value, error_trace = sys.exc_info()
            error_log = "Warning log: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
            errored_partition.append(["{}_{}_{}_{}".format(partition.ref_chrom, partition.ref_start, partition.ref_end, partition.target_variant_type), error_log])
            continue

    interval_list_out.close()
    interval_vcf_out.close()

    interval_end_time = datetime.datetime.now()

    logging.info("Collecting {} {}_{}_{}, {} candidate events found. Time cost {}s".format(options.sample_name, interval_chrom, interval_start, interval_end, len(partitions_list), (interval_end_time - interval_start_time).seconds))

    return errored_interval, errored_partition
