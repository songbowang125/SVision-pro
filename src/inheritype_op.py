import os
import pysam
import torch
from torch.utils.data import DataLoader
from src.network import op_dataset as dataset
from src.network.model_miniunet import MiniUnet
from src.network.model_liteunet import LiteUnet

from src.network.model_unet import Unet
from src.network.utils import image_helpers, image_transforms
from src.output_vcf import output_it_record, classify_vaf_to_genotype
import logging
import sys, traceback
import datetime
import numpy as np


class PIX_CLASS:
    def __init__(self, class_type, pix_start, pix_end):
        self.class_type = class_type
        self.pix_start = pix_start
        self.pix_end = pix_end
        self.occupied_pixes = self.pix_end - self.pix_start + 1

    def to_string(self):
        return "{}-{}-{}".format(self.class_type, self.pix_start, self.pix_end)


class PIX_ALLELE:
    def __init__(self, source, included_classes):

        self.source = source
        self.total_occupied_pixes = sum([clas.occupied_pixes for clas in included_classes])

        self.ref_class = "REF"
        self.ref_vaf = 0
        self.alt_vaf = 0
        self.alt_class = "NA"

        self.included_class_types = {}

        for clas in included_classes:
            cur_class_type = clas.class_type
            cur_class_vaf = round(clas.occupied_pixes / self.total_occupied_pixes, 2)

            if cur_class_type == "REF":
                self.ref_vaf += cur_class_vaf
            else:
                self.alt_vaf += cur_class_vaf
                self.alt_class = cur_class_type


    def to_string(self):

        return "{}: {}-{},{}-{}".format(self.source, self.ref_class, self.ref_vaf, self.alt_class, self.alt_vaf)



def parse_pix_classes_to_alleles(pred_classes):
    """
    for each ref column in pred img, calculate the vaf
    """
    allele_list = []

    # # STEP: for each ref index column, the first (up) and the last (bottom) must bu background
    if pred_classes[0].class_type != "BG" or pred_classes[-1].class_type != "BG":
        return allele_list

    # # STEP: traverse res for the ref column
    next_allele_source = "base"
    included_classes = []
    for i in range(1, len(pred_classes)):

        # # when meeting background, then we find a whole event and output it
        if pred_classes[i].class_type == "BG":
            tmp_allele = PIX_ALLELE(next_allele_source, included_classes)
            # if tmp_vaf.total_occupied_pixes == 17:
            allele_list.append(tmp_allele)

            included_classes = []

            if next_allele_source == "base":
                next_allele_source = "target"
            else:
                next_allele_source = "base"

        else:
            included_classes.append(pred_classes[i])

    return allele_list


def parse_pred_matrix_on_ref_and_read(pred_type_matrix):
    """
    parse predict matrix to classes and freqs
    """

    matrix_size = np.shape(pred_type_matrix)[0]
    alleles_of_ref_indexs = []

    # # STEP: collect predict class for each ref index
    for index_on_ref in range(matrix_size):

        pix_classes = []  # format likes: ['BG-0-59', 'REF-60-84', 'BG-85-86', 'DUP-87-111', 'BG-112-255']

        previous_class = previous_start = previous_end = -1

        # # STEP: traverse each read index to collect corresponding classed by given class colors
        for index_on_read in range(matrix_size):
            # # meet the end of the matrix
            if index_on_read == matrix_size - 1:
                previous_end = index_on_read
                if previous_class != "INS":  # INS was deal on read axis
                    pix_classes.append(PIX_CLASS(previous_class, previous_start, previous_end))
            else:
                # # # find the matched class according to color valudes

                cur_type = pred_type_matrix[index_on_read, index_on_ref]

                if cur_type != previous_class:
                    if previous_class != -1:
                        # append to list
                        if previous_class != "INS":
                            pix_classes.append(PIX_CLASS(previous_class, previous_start, previous_end))

                    # update all info
                    previous_class = cur_type
                    previous_start = index_on_read
                    previous_end = index_on_read

                else:
                    # update end position
                    previous_end = index_on_read

        # # group each position's res to intervals
        pix_alleles = parse_pix_classes_to_alleles(pix_classes)  # # format likes: ['base: REF-1.0,NA-0', 'target: REF-0,DUP-1.0']

        alleles_of_ref_indexs.append(pix_alleles)

    # # ---------------------------------------------------------------------
    # # perform read index
    alleles_of_read_indexs = []

    # # STEP: collect predict class for each ref index
    for index_on_read in range(matrix_size):
        # continue
        pix_classes = []

        previous_class = previous_start = previous_end = -1

        # # STEP: traverse each read index to collect corresponding classed by given class colors
        for index_on_ref in range(matrix_size):
            if index_on_ref == matrix_size - 1:
                previous_end = index_on_ref
                if previous_class in ["INS", "BG", "REF"]:
                    pix_classes.append(PIX_CLASS(previous_class, previous_start, previous_end))
            else:
                # # find the matched class according to color valudes

                cur_type = pred_type_matrix[index_on_read, index_on_ref]

                if cur_type != previous_class:
                    if previous_class != -1:
                        # append to list
                        if previous_class in ["INS", "BG", "REF"]:
                            pix_classes.append(PIX_CLASS(previous_class, previous_start, previous_end))

                    # update all info
                    previous_class = cur_type
                    previous_start = index_on_ref
                    previous_end = index_on_ref

                else:
                    # update end position
                    previous_end = index_on_ref

        # print(index_on_read, [i.to_string() for i in pix_classes], len(pix_classes))

        # # STEP: fix denovo bkp for INS
        if len(pix_classes) == 3 and pix_classes[1].class_type == "INS":
            ins_pix_num = pix_classes[1].pix_end - pix_classes[1].pix_start + 1

            pix_classes[0].pix_end = pix_classes[0].pix_end - ins_pix_num - 1

            tmp_ref = PIX_CLASS("REF", pix_classes[0].pix_end + 1, pix_classes[0].pix_end + ins_pix_num)
            tmp_bg = PIX_CLASS("BG", tmp_ref.pix_end + 1, tmp_ref.pix_end + 2)

            pix_classes.insert(1, tmp_bg)
            pix_classes.insert(1, tmp_ref)

        if len(pix_classes) == 4 and pix_classes[2].class_type == "INS":

            ins_pix_num = pix_classes[2].pix_end - pix_classes[2].pix_start + 1

            pix_classes[0].pix_end = pix_classes[0].pix_end - ins_pix_num - 1

            tmp_ref = PIX_CLASS("REF", pix_classes[0].pix_end + 1, pix_classes[0].pix_end + ins_pix_num)
            tmp_bg = PIX_CLASS("BG", tmp_ref.pix_end + 1, tmp_ref.pix_end + 2)

            pix_classes.insert(1, tmp_bg)
            pix_classes.insert(1, tmp_ref)

        # # group each position's res to intervals
        pix_alleles = parse_pix_classes_to_alleles(pix_classes)

        alleles_of_read_indexs.append(pix_alleles)

    return alleles_of_ref_indexs, alleles_of_read_indexs


def run_unet_predict(data_path, options, sigmod_thresh=0.8):
    """
    run prediction for data
    """

    # start_time = datetime.datetime.now()

    model_path = options.model_path
    device = options.device
    crop_size = options.img_size

    img_alleles_dict = {}

    # # STEP: load data
    #image_norm_transform = image_transforms.ImageNormalize()
    image_norm_transform = None

    center_crop_transform = image_transforms.CenterCrop(crop_size)
    img2tensor_transform = image_transforms.NpyToTensor()
    mask2tensor_transform = image_transforms.MaskToTensor()

    test_set = dataset.Dataset(data_path, 'test', 1, image_norm_transform=image_norm_transform, img2tensor_transform=img2tensor_transform, center_crop_transform=center_crop_transform, mask2tensor_transform=mask2tensor_transform)
    # test_loader = DataLoader(test_set, batch_size=options.batch_size, shuffle=False, num_workers=options.unet_cpu_num)
    test_loader = DataLoader(test_set, batch_size=options.batch_size, shuffle=False)

    palette = dataset.palette
    palette_order = dataset.palette_order
    num_classes = dataset.num_classes

    # # STEP: load model
    if device == "gpu":
        torch.cuda.set_device(options.gpu_id)
        # net = MiniUnet(img_ch=3, num_classes=num_classes, depth=2).cuda()
        # net = Unet(img_ch=3, num_classes=num_classes, depth=2).cuda()
        net = LiteUnet(img_ch=3, num_classes=num_classes, depth=2).cuda()

        net.load_state_dict(torch.load(model_path, map_location=torch.device('cuda')))
        net.eval()
    else:
        torch.set_num_threads(options.unet_cpu_num)
        # net = MiniUnet(img_ch=3, num_classes=num_classes, depth=2)
        #net = Unet(img_ch=3, num_classes=num_classes, depth=2)
        net = LiteUnet(img_ch=3, num_classes=num_classes, depth=2)

        net.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
        net.eval()

    # part_time1 = datetime.datetime.now()
    # print(data_path, (part_time1 - start_time).seconds)

    # # STEP: predict each img
    for input, file_name in test_loader:

        # print(img_name)
        if device == "gpu":
            X = input.cuda()
        else:
            X = input

        # # STEP: generate pred matrix
        pred = net(X)
        #pred = torch.sigmoid(pred)
        pred = pred.cpu().detach()

        #pred[pred < sigmod_thresh] = 0

        for i in range(options.batch_size):
            img_name = file_name[i].split("/")[-1]

            # # STEP: output pred img
            pred_matrix = image_helpers.onehot_to_mask(np.array(pred[i]).squeeze().transpose([1, 2, 0]), palette)
            image_helpers.array_to_img(pred_matrix).save(file_name[0].replace(".png", ".predict.png"))

            # # STEP: parse matrix
            pred_type_matrix = image_helpers.onehot_to_type(np.array(pred[i]).squeeze().transpose([1, 2, 0]), palette_order)

            alleles_of_ref_indexs, alleles_of_read_indexs = parse_pred_matrix_on_ref_and_read(pred_type_matrix)

            img_alleles_dict[img_name] = {"ref_index": alleles_of_ref_indexs, "read_index": alleles_of_read_indexs}

    # part_time2 = datetime.datetime.now()
    # print(data_path, (part_time2 - part_time1).seconds)

    return img_alleles_dict


def classify_all_inheritypes_to_final(simple_var_inheritypes, report_new_bkps=False):
    included_inheritypes = list(set(simple_var_inheritypes))

    # # STEP: only one inheritype for all positions
    if len(included_inheritypes) == 0:
        raise Exception("Skipped")

    elif len(included_inheritypes) == 1:
        return included_inheritypes[0]

    # # STEP: more than one, only determine denovo bkps
    else:

        simple_var_inheritypes_collect = []  # # format: ["somatic", start_index, end_index]
        previous_inheritype = simple_var_inheritypes[0]
        previous_start_index = 0

        for i in range(1, len(simple_var_inheritypes)):

            if simple_var_inheritypes[i] != previous_inheritype:
                simple_var_inheritypes_collect.append([previous_inheritype, previous_start_index, i - 1])

                previous_inheritype = simple_var_inheritypes[i]
                previous_start_index = i

            # # meet the end
            if i == len(simple_var_inheritypes) - 1:
                simple_var_inheritypes_collect.append([previous_inheritype, previous_start_index, i])

        if report_new_bkps:
            # only determine denovo bkps at the start or end
            # if simple_var_inheritypes_collect[-1][0] == "NewCOMP" and simple_var_inheritypes_collect[-2][0] in ["Germline", "NewALLELE"]:
            if simple_var_inheritypes_collect[-1][0] == "NewCOMP" and simple_var_inheritypes_collect[-2][0] in ["Germline"]:
                return "NewBKP"

            # if simple_var_inheritypes_collect[0][0] == "NewCOMP" and simple_var_inheritypes_collect[1][0] in ["Germline", "NewALLELE"]:
            if simple_var_inheritypes_collect[0][0] == "NewCOMP" and simple_var_inheritypes_collect[1][0] in ["Germline"]:
                return "NewBKP"

        else:
            new_comp_ratio = 0
            new_allele_ratio = 0
            germline_ratio = 0

            for inheritype, start_index, end_index in simple_var_inheritypes_collect:
                if inheritype == "NewCOMP":
                    new_comp_ratio += end_index - start_index + 1
                elif inheritype == "NewALLELE":
                    new_allele_ratio += end_index - start_index + 1
                elif inheritype == "Germline":
                    germline_ratio += end_index - start_index + 1
                else:
                    pass
            new_comp_ratio = new_comp_ratio / len(simple_var_inheritypes)
            new_allele_ratio = new_allele_ratio / len(simple_var_inheritypes)
            germline_ratio = germline_ratio / len(simple_var_inheritypes)

            # print(germline_ratio, new_allele_ratio)
            if germline_ratio >= 0.2:
                return "Germline"
            elif new_allele_ratio >= 0.2:
                return "NewALLELE"
            else:
                return "NewCOMP"
        # return collections.Counter(simple_var_inheritypes).most_common()[0][0]


def classify_vaf_to_inheritype(base_vaf, target_vaf, diff_freq, options):
    """
    classify vaf types, including germline, somatic, allele_altered
    """

    if options.detect_mode == "somatic":
        if base_vaf == 0 and target_vaf >= diff_freq:
            return "NewCOMP"
    else:
        if base_vaf == 0:
            return "NewCOMP"

    if abs(base_vaf - target_vaf) >= diff_freq:
        return "NewALLELE"

    return "Germline"


def classify_allele_to_inheritype(base_allele, target_allele, diff_freq, options):
    """
    classify vaf types, including germline, somatic, allele_altered
    """

    if options.detect_mode == "somatic":

        if base_allele.alt_vaf == 0 and target_allele.alt_vaf >= diff_freq:
            return "NewCOMP"

        if base_allele.alt_class != target_allele.alt_class and target_allele.alt_vaf >= diff_freq:
            return "NewCOMP"
    else:
        if base_allele.alt_vaf == 0:
            return "NewCOMP"
        if base_allele.alt_class != target_allele.alt_class:
            return "NewCOMP"

    if base_allele.alt_class == target_allele.alt_class and abs(base_allele.alt_vaf - target_allele.alt_vaf) >= diff_freq:
        return "NewALLELE"

    return "Germline"


def inheritype_variants_in_interval(interval_chrom, interval_start, interval_end, options):
    """
    perform inherityping in each partial interval
    """

    # os.environ["CUDA_VISIBLE_DEVICES"] = options.gpu_id

    interval_start_time = datetime.datetime.now()

    errored_genotypes = []
    errored_interval = []
    passed_genotypes_num = 0

    interval_out_path = os.path.join(options.sample_out_path, "{}_{}_{}".format(interval_chrom, interval_start, interval_end))
    if not os.path.exists(interval_out_path):
        return errored_interval, errored_genotypes

    interval_out_vaf = os.path.join(interval_out_path, "vcf.txt")
    interval_it_out_vaf = open(os.path.join(interval_out_path, "vcf.it.txt"), "w")


    try:
        if not options.skip_inheritype:
            # # there is no img list file
            interval_img_file = os.path.join(interval_out_path, "test.txt")
            if not os.path.exists(interval_img_file):
                return errored_interval, errored_genotypes

            # # the img list file is empty
            interval_img_file = open(interval_img_file)
            interval_img_num = len(interval_img_file.readlines())
            interval_img_file.close()

            if interval_img_num != 0:
                # # STEP: run prediction with trained unet
                img_unet_res_dict = run_unet_predict(interval_out_path, options)
            else:
                img_unet_res_dict = {}

        else:
            img_unet_res_dict = {}
    except:
        error_type, error_value, error_trace = sys.exc_info()
        error_log = "Error log: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
        errored_interval = ["{}_{}_{}".format(interval_chrom, interval_start, interval_end), error_log]

        return errored_interval, errored_genotypes

    # # STEP: match vaf to partition's variant
    for record in pysam.VariantFile(interval_out_vaf):

        try:
            #raise Exception("Image not found")
            left_extension_start = int(record.info["LEFT"])
            left_extension_length = int(record.info["LEFTLEN"])

            color_plot_resize_ratio = float(record.info["RESIZE"])
            simple_var_list = record.info["BKPS"]
            pseudo_vaf = float(record.info["VAF"])
            pseudo_base_vaf = record.info["BASE"]

            # # STEP: assign VAF to target variants
            vars_bkp_list = []

            vars_length_list = []
            vars_inheritype_list = []
            vars_base_gt_list = []
            var_read_start_index = int(left_extension_length / color_plot_resize_ratio)

            for simple_var_index in range(len(simple_var_list)):
                simple_var_split = simple_var_list[simple_var_index].split("++")

                simple_var_type = simple_var_split[0]
                simple_var_hybrid_length = int(simple_var_split[1])
                simple_var_ref_chrom = simple_var_split[2]
                simple_var_ref_start = int(simple_var_split[3])
                simple_var_ref_end = int(simple_var_split[4])
                simple_var_ref_insert = int(simple_var_split[5])

                vars_inheritype_list.append([simple_var_type])
                vars_base_gt_list.append([simple_var_type])
                if options.bkp_mode == "extended":
                    vars_bkp_list.extend([simple_var_ref_start, simple_var_ref_end, simple_var_ref_insert])
                else:
                    if simple_var_type not in ["dDUP", "idDUP"]:
                        vars_bkp_list.extend([simple_var_ref_start, simple_var_ref_end, simple_var_ref_insert])
                    else:
                        vars_bkp_list.extend([simple_var_ref_insert])

                vars_length_list.append(simple_var_hybrid_length)
                var_ref_start_index = int((simple_var_ref_start - left_extension_start) / color_plot_resize_ratio)

                var_len = simple_var_hybrid_length

                var_len_resized = len(np.arange(0, var_len, color_plot_resize_ratio))

                # # avoid list out of range
                if var_ref_start_index < 0:
                    var_ref_start_index = 0
                if var_ref_start_index >= options.img_size or (var_ref_start_index + var_len_resized) > options.img_size:
                    var_ref_start_index = options.img_size - var_len_resized - 1

                if simple_var_type in ["tDUP", "dDUP"]:
                    simple_var_type_tmp = "DUP"
                elif simple_var_type in ["itDUP", "idDUP"]:
                    simple_var_type_tmp = "invDUP"
                else:
                    simple_var_type_tmp = simple_var_type

                # # # this var is out of the boundary of
                # if var_ref_start_index < 0 or var_ref_start_index > len(vaf_of_each_ref_index):
                #     continue
                color_plot_names = record.info["IMG"]

                if not options.skip_inheritype and len(color_plot_names) == 0:
                    raise Exception("Image not found")

                for color_plot_name_index in range(len(color_plot_names)):
                    color_plot_name = color_plot_names[color_plot_name_index]

                    # print(color_plot_name)
                    if color_plot_name in img_unet_res_dict:
                        alleles_of_each_ref_index = img_unet_res_dict[color_plot_name]["ref_index"]
                        alleles_of_each_read_index = img_unet_res_dict[color_plot_name]["read_index"]

                        # previous_target_vaf_val = -1
                        simple_var_inheritypes = []
                        simple_var_target_vaf = []
                        simple_var_base_vaf = []

                        candidate_pos = 0

                        for i in np.arange(0, var_len, color_plot_resize_ratio):
                            # print(len(alleles_of_each_read_index), i ,candidate_pos, var_read_start_index, var_read_start_index + candidate_pos)
                            if simple_var_type in ["INS"]:
                                base_target_allele = alleles_of_each_read_index[var_read_start_index + candidate_pos]
                            else:
                                base_target_allele = alleles_of_each_ref_index[var_ref_start_index + candidate_pos]

                            candidate_pos += 1

                            if len(base_target_allele) != 0 and len(base_target_allele) % 2 == 0 and base_target_allele[1].alt_class == simple_var_type_tmp:
                                base_allele = base_target_allele[0]
                                target_allele = base_target_allele[1]

                                simple_var_inheritypes.append(classify_allele_to_inheritype(base_allele, target_allele, options.min_diff_freq, options))

                                simple_var_target_vaf.append(target_allele.alt_vaf)
                                simple_var_base_vaf.append(base_allele.alt_vaf)

                        simple_var_final_inheritype = classify_all_inheritypes_to_final(simple_var_inheritypes, options.report_new_bkps)

                        # # STEP: calculate the base GT
                        if simple_var_final_inheritype == "NewCOMP":
                            # final_base_vaf = float(pseudo_base_vaf[color_plot_name_index])
                            final_base_vaf = 0.0

                        elif simple_var_final_inheritype == "NewALLELE":
                            final_base_vaf = []
                            for i in range(len(simple_var_inheritypes)):
                                if simple_var_inheritypes[i] == "NewALLELE":
                                    final_base_vaf.append(simple_var_base_vaf[i])

                            final_base_vaf = np.average(final_base_vaf)
                        elif simple_var_final_inheritype == "Germline":
                            final_base_vaf = []
                            for i in range(len(simple_var_inheritypes)):
                                if simple_var_inheritypes[i] == "Germline":
                                    final_base_vaf.append(simple_var_base_vaf[i])

                            final_base_vaf = np.average(final_base_vaf)
                        else:
                            final_base_vaf = 0.0

                        if simple_var_final_inheritype in ["Germline", "NewALLELE", "NewBKP"] and final_base_vaf != float(pseudo_base_vaf[color_plot_name_index]):
                            final_base_vaf = float(pseudo_base_vaf[color_plot_name_index])

                        if simple_var_final_inheritype in ["NewCOMP"] and final_base_vaf != float(pseudo_base_vaf[color_plot_name_index]):
                            raise Exception("Inconsistent inheritype")

                        base_gt = classify_vaf_to_genotype(final_base_vaf)

                        vars_base_gt_list[simple_var_index].append(base_gt)
                        vars_inheritype_list[simple_var_index].append(simple_var_final_inheritype)

                    else:
                        raise Exception("Image not found")

                if simple_var_type != "DEL":
                    var_read_start_index += var_len_resized

                vars_inheritype_list[simple_var_index] = "_".join(vars_inheritype_list[simple_var_index])
            output_it_record(record, vars_length_list, vars_bkp_list, vars_base_gt_list, vars_inheritype_list, interval_it_out_vaf, pseudo_vaf, options.bkp_mode)
            passed_genotypes_num += 1

        except:

            # # STEP: assign VAF to target variants
            vars_bkp_list = []
            vars_length_list = []
            vars_base_gt_list = []
            vars_inheritype_list = []

            simple_var_list = record.info["BKPS"]
            pseudo_vaf = float(record.info["VAF"])

            if len(options.base_path) != 0:
                pseudo_base_gt = [gt_info.split(":")[0] for gt_info in str(record).strip().split("\t")[10: ]]
                pseudo_base_vaf = record.info["BASE"]
            else:
                pseudo_base_gt = []
                pseudo_base_vaf = []

            for simple_var_index in range(len(simple_var_list)):
                simple_var_split = simple_var_list[simple_var_index].split("++")

                simple_var_type = simple_var_split[0]
                simple_var_hybrid_length = int(simple_var_split[1])
                simple_var_ref_chrom = simple_var_split[2]
                simple_var_ref_start = int(simple_var_split[3])
                simple_var_ref_end = int(simple_var_split[4])
                simple_var_ref_insert = int(simple_var_split[5])

                vars_length_list.append(simple_var_hybrid_length)

                vars_base_gt_list.append([simple_var_type])
                vars_inheritype_list.append([simple_var_type])

                if len(options.base_path) != 0:
                    vars_base_gt_list[simple_var_index].extend(pseudo_base_gt)
                    for i in range(len(options.base_path)):
                        vars_inheritype_list[simple_var_index].append(classify_vaf_to_inheritype(float(pseudo_base_vaf[i]), pseudo_vaf, options.min_diff_freq, options))

                if options.bkp_mode == "extended":
                    vars_bkp_list.extend([simple_var_ref_start, simple_var_ref_end, simple_var_ref_insert])
                else:
                    if simple_var_type not in ["dDUP", "idDUP"]:
                        vars_bkp_list.extend([simple_var_ref_start, simple_var_ref_end, simple_var_ref_insert])
                    else:
                        vars_bkp_list.extend([simple_var_ref_insert])

                vars_inheritype_list[simple_var_index] = "_".join(vars_inheritype_list[simple_var_index])

            output_it_record(record, vars_length_list, vars_bkp_list, vars_base_gt_list, vars_inheritype_list, interval_it_out_vaf, pseudo_vaf, options.bkp_mode)

            error_type, error_value, error_trace = sys.exc_info()
            error_log = "Warning log: " + str(error_value) + '. ' + 'Locate At: ' + str(traceback.extract_tb(error_trace))
            errored_genotypes.append(["{}_{}_{}_{}".format(record.contig, record.start, record.stop, record.info["SVTYPE"]), error_log])

            continue

    interval_end_time = datetime.datetime.now()
    if not options.skip_inheritype:
        logging.info("Predicting {} {}_{}_{}, {} events pass lite-Unet module. Time cost: {}s".format(options.sample_name, interval_chrom, interval_start, interval_end, passed_genotypes_num + len(errored_genotypes), (interval_end_time - interval_start_time).seconds))

    interval_it_out_vaf.close()

    return errored_interval, errored_genotypes

