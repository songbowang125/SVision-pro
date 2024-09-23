import numpy as np
from src.version import __version__
import pysam
import os
import logging
import datetime


def binomial_probability(k,n,p):
    try:
        #Binomial coef cancels out for likelihood ratios
        #return binomial_coef(n,k) * (p**k) * ((1.0-p)**(n-k))
        return (p**k) * ((1.0-p)**(n-k))
    except OverflowError:
        return 1.0


def classify_supp_to_genotype(support, coverage):
    genotype_p = [((0, 0), 0.05),
                  ((0, 1), 0.5),
                  ((1, 1), 0.95)]

    normalized_support = 16
    normalized_coverage = 20

    genotype_likelihoods = []
    for gt, p in genotype_p:
        q = binomial_probability(support, coverage, p)
        genotype_likelihoods.append((gt, q))
    genotype_likelihoods.sort(key=lambda k: k[1], reverse=True)

    sum_likelihoods = sum(q for gt, q in genotype_likelihoods)
    normalized_likelihoods = [(gt, (q / sum_likelihoods)) for gt, q in genotype_likelihoods]

    gt1, q1 = normalized_likelihoods[0]
    gt2, q2 = normalized_likelihoods[1]

    print(gt1)


def classify_vaf_to_genotype(vaf):
    """
    classify genotype to 1/0, 0/1, 0/0
    """

    if vaf >= 0.8:
        return '1/1'
    elif vaf >= 0.2:
        return '0/1'
    else:
        return '0/0'

def re_merge_bnds(raw_bnd_records, options):
    """
    due the primary and supplementary, the bnd start could be either in the first chrom or the second chrom.
    In that case, we need merge the two condition together
    """

    def parse_bnd_record(bnd_record, reverse_two_bkp=False):

        reverse_orient = {"+": "-", "-": "+"}

        first_chrom, first_pos, first_orient = bnd_record.contig, bnd_record.start + 1, "+"
        second_info = list(bnd_record.alts)[0]

        if "]" in second_info:
            second_orient = "-"
        else:
            second_orient = "+"

        second_info = second_info.replace("[N", "").replace("]N", "").replace("]", "").replace("[", "").split(":")
        second_chrom, second_pos = second_info[0], int(second_info[1])

        if reverse_two_bkp:
            return second_chrom, second_pos, reverse_orient[second_orient], first_chrom, first_pos, reverse_orient[first_orient]
        else:
            return first_chrom, first_pos, first_orient, second_chrom, second_pos, second_orient

    # # traverse each
    merged_bnd_records = []

    for base_record in raw_bnd_records:

        merge_flag = False
        for target_record in merged_bnd_records:

            base_first_chrom, base_first_pos, base_first_orient, base_second_chrom, base_second_pos, base_second_orient = parse_bnd_record(base_record)

            target_first_chrom, target_first_pos, target_first_orient, target_second_chrom, target_second_pos, target_second_orient = parse_bnd_record(target_record, reverse_two_bkp=True)

            if base_first_chrom == target_first_chrom and base_second_chrom == target_second_chrom and base_first_orient == target_first_orient and base_second_orient == target_second_orient and abs(base_first_pos - target_first_pos) < options.dist_diff_length and abs(base_second_pos - target_second_pos) < options.dist_diff_length:
                merge_flag = True
                break

        if not merge_flag:
            merged_bnd_records.append(base_record)

    return merged_bnd_records



def output_vcf_final(output_vcf_path, intervals_list, options):
    """
    merge split gt-vcf files to the final one
    """

    start_time = datetime.datetime.now()

    output_cnt = 0
    csv_cnt = 0
    sv_cnt = 0

    with open(output_vcf_path, "w") as fout:

        output_vcf_header_final(options.genome_path, options.target_path, options.base_path, fout)

        raw_bnds_records = []

        for interval_chrom, interval_start, interval_end in intervals_list:
            interval_info = "{}_{}_{}".format(interval_chrom, interval_start, interval_end)

            interval_it_vcf = os.path.join(options.sample_out_path, "{}/vcf.it.txt".format(interval_info))

            if not os.path.exists(interval_it_vcf):
                continue

            with pysam.VariantFile(interval_it_vcf) as interval_records:
                for record in interval_records:

                    if not options.skip_bnd and record.info["SVTYPE"] == "BND":
                        raw_bnds_records.append(record)
                        continue

                    record_split = str(record).split("\t")

                    # # update the id
                    record_split[2] = str(output_cnt)
                    output_cnt += 1

                    if record_split[4] == "SV":
                        sv_cnt += 1
                    else:
                        csv_cnt += 1

                    record_split[7] = record_split[7].replace("++", "_")

                    fout.write("\t".join(record_split))

        if not options.skip_bnd:
            merged_bnd_records = re_merge_bnds(raw_bnds_records, options)

            for record in merged_bnd_records:
                record.info.pop("SVLEN")
                record.info.pop("BKPS")

                fout.write(str(record))

    end_time = datetime.datetime.now()
    logging.info("Outputting {} to {}. Time cost: {}s".format(options.sample_name, output_vcf_path, (end_time - start_time).seconds))

    return output_cnt, sv_cnt, csv_cnt


def output_vcf_header_final(ref_path, target_path, base_path, output_vcf):
    """
    generate vcf header
    """

    target_name = target_path.split("/")[-1]

    base_names = []
    for base_name in base_path:
        base_names.append(base_name.split("/")[-1])
    base_names = "\t".join(base_names)

    print("##fileformat=VCFv4.3", file=output_vcf)
    print("##source=SVision v{0}".format(__version__), file=output_vcf)

    # # STEP: add chromosome info
    ref_file = pysam.FastaFile(ref_path)
    chroms = ref_file.references
    for chr in chroms:
        chr_length = ref_file.get_reference_length(chr)
        print('##contig=<ID={0},length={1}>'.format(chr, chr_length), file=output_vcf)

    # # STEP: add info field
    print("##CHROM=<CHROM=XXX,Description=\"Chromosome ID\">", file=output_vcf)
    print("##POS=<POS=XXX,Description=\"Start position of the SV described in this region\">", file=output_vcf)
    print("##ID=<ID=XXX,Description=\"ID of the SV described in this region\">", file=output_vcf)
    print("##REF=<REF=N,Description=\"Ref's sequence in that region, default=N\">", file=output_vcf)
    print("##QUAL=<QUAL=XXX,Description=\"The SV quality of the SV described in this region\">", file=output_vcf)

    print("##ALT=<ID=SV,Description=\"Simple SVs\">", file=output_vcf)
    print("##ALT=<ID=CSV,Description=\"Complex or nested SVs\">", file=output_vcf)

    print("##FILTER=<ID=PASS,Description=\"Passed\">", file=output_vcf)
    print("##FILTER=<ID=POLY,Description=\"Polymorphic\">", file=output_vcf)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV described in this region\">", file=output_vcf)
    print("##INFO=<ID=IT,Number=1,Type=String,Description=\"Inheritype of the whole event \">", file=output_vcf)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=output_vcf)
    print("##INFO=<ID=BKPS,Number=.,Type=String,Description=\"All breakpoints (length-start-end-insert) in this region.\">", file=output_vcf)
    print("##INFO=<ID=BKPSIT,Number=.,Type=String,Description=\"All breakpoints' inheritype in this region.\">", file=output_vcf)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"CNN predicted SV type, containing INS, DEL, DUP, tDUP (tandem duplication) and INV\">", file=output_vcf)
    print("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"SV support number in this region\">", file=output_vcf)
    print("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"SV allele frequency in this region\">", file=output_vcf)
    print("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"SV support read names in this region\">", file=output_vcf)

    # # STEP: add gt info
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=output_vcf)
    print("##FORMAT=<ID=DR,Number=1,Type=String,Description=\"high-quality reference reads\">", file=output_vcf)
    print("##FORMAT=<ID=DV,Number=1,Type=String,Description=\"high-quality variant reads\">", file=output_vcf)

    if len(base_path) == 0:
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{target_name}", file=output_vcf)
    else:
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{target_name}\t{base_names}", file=output_vcf)


def output_it_record(origin_record, vars_length_list, var_bkp_list, vars_base_gt_list, vars_inheritype_list, output_vcf, pseudo_vaf, bkp_mode="affected"):
    """
    output ingerityped record to vcf file
    """

    origin_record_str_split = str(origin_record).strip().split("\t")
    origin_gt_info_list = origin_record_str_split[9: ]

    final_ref_start = min(var_bkp_list)
    final_ref_end = max(var_bkp_list)

    if bkp_mode == "affected":
        sv_len = sum(vars_length_list)
    else:
        if origin_record.info["SVTYPE"] == "INS":
            sv_len = sum(vars_length_list)
        else:
            sv_len = final_ref_end - final_ref_start + 1

    target_gt = classify_vaf_to_genotype(pseudo_vaf)

    base_file_num = len(vars_base_gt_list[0]) - 1

    target_dr = origin_gt_info_list[0].split(":")[1]
    target_dv = origin_gt_info_list[0].split(":")[2]

    if base_file_num == 0:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};RNAMES={}\tGT\t{}".format(origin_record.contig, final_ref_start, 0, "N", origin_record.alts[0], ".", origin_record.filter.keys()[0], final_ref_end, sv_len, origin_record.info["SVTYPE"], origin_record.info["SUPPORT"], pseudo_vaf, ",".join(origin_record.info["BKPS"]), ",".join(origin_record.info["RNAMES"]), target_gt), file=output_vcf)
    else:
        #base_origin_gt_info = origin_record_str_split[-base_file_num:]

        base_gt_list = []
        base_dr_list = []
        base_dv_list = []

        for base_file_index in range(1, len(vars_base_gt_list[0])):
            base_gt_list.append([])
            for simple_var_index in range(len(vars_base_gt_list)):
                base_gt_list[base_file_index - 1].append(vars_base_gt_list[simple_var_index][base_file_index])

            base_dr_list.append(origin_gt_info_list[base_file_index].split(":")[1])
            base_dv_list.append(origin_gt_info_list[base_file_index].split(":")[2])

        # # base gt preference
        for base_file_index in range(len(base_gt_list)):
            if "1/1" in base_gt_list[base_file_index]:
                base_gt_list[base_file_index] = "1/1"
            elif "0/1" in base_gt_list[base_file_index]:
                base_gt_list[base_file_index] = "0/1"
            elif "0/0" in base_gt_list[base_file_index]:
                base_gt_list[base_file_index] = "0/0"
            else:
                base_gt_list[base_file_index] = "./."

            # # refix gt by inheritype
            #if "Germline" in ",".join(vars_inheritype_list):
            #    base_gt_list[base_file_index] = target_gt

            base_file_inheritypes = [var_it.split("_")[base_file_index + 1] for var_it in vars_inheritype_list]
            if "Germline" in base_file_inheritypes:
                base_gt_list[base_file_index] = target_gt

        base_gt = "\t".join(["{}:{}:{}".format(base_gt_list[i], base_dr_list[i], base_dv_list[i]) for i in range(len(base_gt_list))])


        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};BKPSIT={};RNAMES={}\tGT:DR:DV\t{}:{}:{}\t{}".format(origin_record.contig, final_ref_start, 0, "N", origin_record.alts[0], ".", origin_record.filter.keys()[0], final_ref_end, sv_len, origin_record.info["SVTYPE"], origin_record.info["SUPPORT"], pseudo_vaf, ",".join(origin_record.info["BKPS"]), ",".join(vars_inheritype_list), ",".join(origin_record.info["RNAMES"]), target_gt, target_dr, target_dv, base_gt), file=output_vcf)


def determine_pseudo_base_vaf_and_gt(base_partition_dict):
    """
    calculate pseudo base gt if base files are provided.
    For each base file (base_index), we use the first base partition to determine pseudo base gt
    """

    base_vaf_list = []
    base_gt_list = []
    base_dv_list = []
    base_dr_list = []

    for base_index in base_partition_dict:

        if len(base_partition_dict[base_index]) == 0:
            base_vaf = 0.0
            base_dv = "NA"
            base_dr = "NA"
        else:
            base_vaf = base_partition_dict[base_index][0].included_hyper_cigars[0].calculate_pseudo_vaf()

            base_dv = len(base_partition_dict[base_index][0].supp_reads)
            base_dr = max(base_partition_dict[base_index][0].included_hyper_cigars[0].calculate_coverage() - base_dv, 0)

        base_vaf_list.append(str(base_vaf))
        base_gt_list.append(classify_vaf_to_genotype(base_vaf))

        base_dv_list.append(base_dv)
        base_dr_list.append(base_dr)

    return base_vaf_list, base_gt_list, base_dv_list, base_dr_list



def output_origin_record_to_vcf(partition, output_vcf, options):
    """
    generate vcf record for each partition in this interval
    """

    supp_reads_str = ",".join(partition.supp_reads)
    supp_reads_num = len(partition.supp_reads)

    # # STEP: define sv and csv
    pseudo_sv_type = partition.target_variant_type

    pseudo_vaf = partition.included_hyper_cigars[0].calculate_pseudo_vaf()

    pseudo_gt = classify_vaf_to_genotype(pseudo_vaf)

    pseudo_dv = supp_reads_num
    pseudo_dr = max(partition.included_hyper_cigars[0].calculate_coverage() - pseudo_dv, 0)

    if len(options.base_path) != 0:
        base_pseudo_vaf_list, base_pseudo_gt_list, base_pseudo_dv_list, base_pseudo_dr_list = determine_pseudo_base_vaf_and_gt(partition.base_partitions)
    else:
        base_pseudo_vaf_list = []
        base_pseudo_gt_list = []
        base_pseudo_dv_list = []
        base_pseudo_dr_list = []

    if "+" in pseudo_sv_type:
        is_csv = "CSV"
    elif pseudo_sv_type in ["itDUP", "idDUP"]:
        is_csv = "CSV"
    else:
        is_csv = "SV"

    # # STEP: define BKP
    bkp_str = []
    for simple_var in partition.target_variant:
        if simple_var.type in ["dDUP", "tDUP", "itDUP", "idDUP"]:
            bkp_str.append("{}++{}++{}++{}++{}++{}".format(simple_var.type, simple_var.hybrid_length, simple_var.dup_source_chrom, simple_var.dup_source_start, simple_var.dup_source_end, simple_var.ref_start))
        else:
            bkp_str.append("{}++{}++{}++{}++{}++{}".format(simple_var.type, simple_var.hybrid_length, simple_var.ref_chrom, simple_var.ref_start, simple_var.ref_end, simple_var.ref_start))

    bkp_str = ",".join(bkp_str)

    if not options.skip_bnd and pseudo_sv_type == "BND":
        simple_var = partition.target_variant[0]

        if simple_var.bnd_source_strand == "+":
            bnd_str = "[{}:{}[N".format(simple_var.bnd_source_chrom, simple_var.bnd_source_start)
        else:
            bnd_str = "]{}:{}]N".format(simple_var.bnd_source_chrom, simple_var.bnd_source_start)

        if len(options.base_path) != 0:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};RNAMES={};LEFT={};LEFTLEN={};IMG={};RESIZE={};BASE={}\tGT:DR:DV\t{}:{}:{}\t{}".format(partition.ref_chrom, partition.ref_start, 0, "N", bnd_str, "NA", partition.filter_marker, partition.ref_end, partition.hybrid_length, pseudo_sv_type, supp_reads_num, partition.included_hyper_cigars[0].calculate_pseudo_vaf(), bkp_str, supp_reads_str, partition.left_extension_start, partition.ref_start - partition.left_extension_start, partition.color_plot_name, partition.color_plot_resize_ratio, ",".join(base_pseudo_vaf_list), pseudo_gt, pseudo_dr, pseudo_dv, "\t".join(["{}:{}:{}".format(base_pseudo_gt_list[i], base_pseudo_dr_list[i], base_pseudo_dv_list[i]) for i in range(len(base_pseudo_gt_list))])), file=output_vcf)
        else:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};RNAMES={};LEFT={};LEFTLEN={};IMG={};RESIZE={}\tGT:DR:DV\t{}:{}:{}".format(partition.ref_chrom, partition.ref_start, 0, "N", bnd_str, "NA", partition.filter_marker, partition.ref_end, partition.hybrid_length, pseudo_sv_type, supp_reads_num, partition.included_hyper_cigars[0].calculate_pseudo_vaf(), bkp_str, supp_reads_str, partition.left_extension_start, partition.ref_start - partition.left_extension_start, partition.color_plot_name, partition.color_plot_resize_ratio, pseudo_gt, pseudo_dr, pseudo_dv), file=output_vcf)

    else:
        # # STEP: output to file
        if len(options.base_path) != 0:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};RNAMES={};LEFT={};LEFTLEN={};IMG={};RESIZE={};BASE={}\tGT:DR:DV\t{}:{}:{}\t{}".format(partition.ref_chrom, partition.ref_start, 0, "N", is_csv, "NA", partition.filter_marker, partition.ref_end, partition.hybrid_length, pseudo_sv_type, supp_reads_num, partition.included_hyper_cigars[0].calculate_pseudo_vaf(), bkp_str, supp_reads_str, partition.left_extension_start, partition.ref_start - partition.left_extension_start, ",".join(partition.color_plot_name), partition.color_plot_resize_ratio, ",".join(base_pseudo_vaf_list), pseudo_gt, pseudo_dr, pseudo_dv, "\t".join(["{}:{}:{}".format(base_pseudo_gt_list[i], base_pseudo_dr_list[i], base_pseudo_dv_list[i]) for i in range(len(base_pseudo_gt_list))])), file=output_vcf)
        else:
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tEND={};SVLEN={};SVTYPE={};SUPPORT={};VAF={};BKPS={};RNAMES={};LEFT={};LEFTLEN={};IMG={};RESIZE={}\tGT:DR:DV\t{}:{}:{}".format(partition.ref_chrom, partition.ref_start, 0, "N", is_csv, "NA", partition.filter_marker, partition.ref_end, partition.hybrid_length, pseudo_sv_type, supp_reads_num, partition.included_hyper_cigars[0].calculate_pseudo_vaf(), bkp_str, supp_reads_str, partition.left_extension_start, partition.ref_start - partition.left_extension_start, ",".join(partition.color_plot_name), partition.color_plot_resize_ratio, pseudo_gt, pseudo_dr, pseudo_dv), file=output_vcf)


def output_vcf_header_in_interval(ref_path, target_path, base_path, output_vcf):
    """
    generate vcf header
    """
    target_name = target_path.split("/")[-1]

    base_names = []
    for base_name in base_path:
        base_names.append(base_name.split("/")[-1])
    base_names = "\t".join(base_names)

    print("##fileformat=VCFv4.3", file=output_vcf)
    print("##source=SVision v{0}".format(__version__), file=output_vcf)

    # # STEP: add chromosome info
    ref_file = pysam.FastaFile(ref_path)
    chroms = ref_file.references
    for chr in chroms:
        chr_length = ref_file.get_reference_length(chr)
        print('##contig=<ID={0},length={1}>'.format(chr, chr_length), file=output_vcf)

    # # STEP: add info field
    print("##CHROM=<CHROM=XXX,Description=\"Chromosome ID\">", file=output_vcf)
    print("##POS=<POS=XXX,Description=\"Start position of the SV described in this region\">", file=output_vcf)
    print("##ID=<ID=XXX,Description=\"ID of the SV described in this region\">", file=output_vcf)
    print("##REF=<REF=N,Description=\"Ref's sequence in that region, default=N\">", file=output_vcf)
    print("##QUAL=<QUAL=XXX,Description=\"The SV quality of the SV described in this region\">", file=output_vcf)

    print("##ALT=<ID=SV,Description=\"Simple SVs\">", file=output_vcf)
    print("##ALT=<ID=CSV,Description=\"Complex or nested SVs\">", file=output_vcf)

    print("##FILTER=<ID=PASS,Description=\"Passed\">", file=output_vcf)
    print("##FILTER=<ID=POLY,Description=\"Polymorphic\">", file=output_vcf)

    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV described in this region\">", file=output_vcf)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=output_vcf)
    print("##INFO=<ID=BKPS,Number=.,Type=String,Description=\"All breakpoints (length-start-end) in this region, where CSV might contain multiple breakpoints.\">", file=output_vcf)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"CNN predicted SV type, containing INS, DEL, DUP, tDUP (tandem duplication) and INV\">", file=output_vcf)
    print("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"SV support number in this region\">", file=output_vcf)
    print("##INFO=<ID=VAF,Number=1,Type=String,Description=\"SV allele frequency in this region\">", file=output_vcf)
    print("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"SV support read names in this region\">", file=output_vcf)

    print("##INFO=<ID=LEFT,Number=1,Type=String,Description=\"Plot parameter, the left extension start\">", file=output_vcf)
    print("##INFO=<ID=LEFTLEN,Number=1,Type=String,Description=\"Plot parameter, the left extension length\">", file=output_vcf)
    print("##INFO=<ID=IMG,Number=.,Type=String,Description=\"Plot parameter, the output image name\">", file=output_vcf)
    print("##INFO=<ID=RESIZE,Number=1,Type=String,Description=\"Plot parameter, the resize ratio\">", file=output_vcf)
    print("##INFO=<ID=BASE,Number=.,Type=String,Description=\"Pseudo base VAF.\">", file=output_vcf)

    # # STEP: add gt info
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=output_vcf)
    print("##FORMAT=<ID=DR,Number=1,Type=String,Description=\"high-quality reference reads\">", file=output_vcf)
    print("##FORMAT=<ID=DV,Number=1,Type=String,Description=\"high-quality variant reads\">", file=output_vcf)

    if len(base_path) == 0:
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{target_name}", file=output_vcf)
    else:
        print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{target_name}\t{base_names}", file=output_vcf)



