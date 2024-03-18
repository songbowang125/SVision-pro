import os
import cv2 as cv
import numpy as np

op_order = ["A", "T", "C", "G", "I", "D", "V"]

# OP_COLOR_DICT = {"M": [0, 0, 0],
#                  "A": [135, 206, 255],
#                   "T": [135, 238, 238],
#                   "C": [135, 186, 150],
#                   "G": [135, 155, 155],
#                  "N": [227, 227, 227],
#                  "REF": [66, 146, 197],
#                  "I": [251, 246, 180],
#                  "D": [253, 229, 217],
#                  "DP": [235, 150, 69],
#                  "V": [127, 206, 187],
#                  }
OP_COLOR_DICT = {"M": [0, 0, 0],
                 "A": [135, 206, 255],
                  "T": [135, 206, 255],
                  "C": [135, 206, 255],
                  "G": [135, 206, 255],
                 "N": [227, 227, 227],
                 "REF": [66, 146, 197],
                 "I": [251, 246, 180],
                 "D": [253, 229, 217],
                 "DP": [235, 150, 69],
                 "V": [127, 206, 187],
                 }
# # convert BGR to RGB

for op in OP_COLOR_DICT:

    tmp = OP_COLOR_DICT[op][0]
    OP_COLOR_DICT[op][0] = OP_COLOR_DICT[op][2]
    OP_COLOR_DICT[op][2] = tmp


match_channel = 2
inv_channel = 1
dup_channel = 0


line_width = 1
gap_width = 1
bar_width = 25

class Line:
    def __init__(self, ref_start, ref_end, read_start, read_end, strand):
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.strand = strand

    def to_string(self):
        return "{}_{}_{}_{}_{}".format(self.ref_start, self.ref_end, self.read_start, self.read_end, self.strand)


def draw_single_line(target_img, line_ref_start, line_ref_end, line_read_start, line_read_end, line_strand):
    """
    draw a single line
    """

    # line_ref_start = int(line_ref_start / resize_ratio)
    # line_ref_end = int(line_ref_end / resize_ratio)
    # line_read_start = int(line_read_start / resize_ratio)
    # line_read_end = int(line_read_end / resize_ratio)

    # # # in case that the line is out of the scope of target img
    # target_img_length_on_ref_axis = np.shape(target_img)[1]
    # line_ref_span = line_ref_end - line_ref_start
    # if line_ref_start > target_img_length_on_ref_axis:
    #     line_ref_end = target_img_length_on_ref_axis
    #     line_ref_start = line_ref_end - line_ref_span

    # # draw lines with different strands
    if line_strand == "+":
        cv.line(target_img, (line_ref_start, line_read_start),
                            (line_ref_end, line_read_end), OP_COLOR_DICT["M"], line_width)
    else:
        cv.line(target_img, (line_ref_start, line_read_end),
                            (line_ref_end, line_read_start), OP_COLOR_DICT["M"], line_width)


def draw_single_cigar_map(target_img, drawable_map, ref_pointer, read_pointer, direction):
    """
    draw a single cigar map at one position
    :return:
    """
    # ref_pointer = int(ref_pointer / resize_ratio)
    # read_pointer = int(read_pointer / read_pointer)

    # # draw base cigar map this map
    if direction == "up":
        read_pointer_going = read_pointer

        for op, channel, freq in drawable_map:

            # op_bar_width = bar_width * (freq / 100)
            # if op_bar_width < 1:
            #     continue
            # else:
            #     op_bar_width = round(op_bar_width)
            op_bar_width = int(round(bar_width * (freq / 100)))

            target_img[max(0, read_pointer_going - op_bar_width): read_pointer_going, ref_pointer, :] = get_drawable_op_color(op, channel)   # # why max? because after subtraction, it will exceed the img size

            read_pointer_going -= op_bar_width

    # # draw target cigar map this map
    elif direction == "down":
        read_pointer_going = read_pointer

        for op, channel, freq in drawable_map:

            op_bar_width = int(round(bar_width * (freq / 100)))

            target_img[read_pointer_going: read_pointer_going + op_bar_width, ref_pointer, :] = get_drawable_op_color(op, channel)
            read_pointer_going += op_bar_width

    elif direction == "left":
        ref_pointer_going = ref_pointer

        for op, channel, freq in drawable_map:
            op_bar_width = int(round(bar_width * (freq / 100)))

            target_img[read_pointer, ref_pointer_going - op_bar_width: ref_pointer_going, :] = get_drawable_op_color(op, channel)
            ref_pointer_going -= op_bar_width

    elif direction == "right":
        ref_pointer_going = ref_pointer

        for op, channel, freq in drawable_map:
            op_bar_width = int(round(bar_width * (freq / 100)))

            target_img[read_pointer, ref_pointer_going: ref_pointer_going + op_bar_width, :] = get_drawable_op_color(op, channel)
            ref_pointer_going += op_bar_width

    else:
        pass


def get_drawable_op_color(op, channel):
    """
    generate drawable op color with different channels
    """

    op = op.upper()

    op_color = OP_COLOR_DICT[op].copy()

    if channel == "inv_channel":
        op_color[inv_channel] -= 100

    elif channel == "dup_channel":
        op_color[dup_channel] -= 100

    elif channel == "invdup_channel":
        op_color[inv_channel] -= 100
        op_color[dup_channel] -= 100

    # # change channel values that are smaller than 0 to 0
    for i in range(len(op_color)):
        if op_color[i] < 0:
            op_color[i] = 0

    return op_color


def get_drawable_map(cigar_map, ref_cigar_map, spec_cigar_map_op, inserted_base=None, mappedI=None):
    """
    get drawable map
    """

    drawable_map = []

    # # there is no specific cigar_map_op
    if spec_cigar_map_op is None:
        spec_cigar_map_op = cigar_map.get_most_frequent_op()

    ref_map_op = ref_cigar_map.get_most_frequent_op()

    # # STEP: calculate rest ref freq
    map_op_freq = cigar_map.freq_map[spec_cigar_map_op]
    rest_ref_freq = 100 - map_op_freq

    if rest_ref_freq > 0:
        drawable_map.append([ref_map_op, "match_channel", rest_ref_freq])

    if spec_cigar_map_op == "D":
        drawable_map.append(["N", "match_channel", map_op_freq])

    elif spec_cigar_map_op == "V":
        drawable_map.append([ref_map_op, "inv_channel", map_op_freq])

    elif spec_cigar_map_op == "I":
        if mappedI is None:
            # # Reset rest ref as N
            if rest_ref_freq > 0:
                drawable_map[0][0] = "N"

            if inserted_base is None:
                drawable_map.append([ref_map_op, "match_channel", map_op_freq])
                # drawable_map.append(["N", "match_channel", map_op_freq])
            else:
                drawable_map.append([inserted_base, "match_channel", map_op_freq])
                # drawable_map.append(["N", "match_channel", map_op_freq])

        elif mappedI == "mapped":
            drawable_map.append([ref_map_op, "dup_channel", map_op_freq])
        elif mappedI == "inv_mapped":
            drawable_map.append([ref_map_op, "invdup_channel", map_op_freq])
        else:
            pass

    elif spec_cigar_map_op in ["A", "T", "C", "G"]:
        drawable_map.append([ref_map_op, "match_channel", map_op_freq])

    return drawable_map


def draw_lines_and_cigar_maps(target_img, left_extension_start, right_extension_start, left_extension_length, right_extension_length, detail_cigars, target_cigar_maps, base_cigar_maps, ref_cigar_maps, resize_ratio):
    """
    draw lines and cigar maps
    """

    # # draw parameters
    img_length_on_read_axis, img_length_on_ref_axis, _ = np.shape(target_img)

    # # STEP: set a global read pointer
    read_pointer_for_draw = 0
    lines = []

    # # STEP: draw part 1, draw left extensions, including both line and maps
    ref_pointer_for_draw = 0
    line_length = 0

    for i in np.arange(0, left_extension_length, resize_ratio):
        target_img[max(read_pointer_for_draw - gap_width - bar_width, 0): read_pointer_for_draw - gap_width, ref_pointer_for_draw, :] = OP_COLOR_DICT["N"]
        target_img[read_pointer_for_draw + gap_width: read_pointer_for_draw + gap_width + bar_width, ref_pointer_for_draw, :] = OP_COLOR_DICT["N"]

        ref_pointer_for_draw += 1
        read_pointer_for_draw += 1

        line_length += 1

    lines.append([ref_pointer_for_draw - line_length, ref_pointer_for_draw, read_pointer_for_draw - line_length, read_pointer_for_draw, "+"])

    # # STEP: draw part 2, draw cigars, including both line and maps
    for detail_cigar in detail_cigars:

        # # STEP: deal with different cigar op
        if detail_cigar.op == "D":

            ref_pointer_for_map = detail_cigar.ref_start - left_extension_start
            ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

            # # set cigar's plot info
            detail_cigar.set_plot_info(ref_pointer_for_draw, read_pointer_for_draw, resize_ratio)

            # # draw map
            for i in range(len(np.arange(0, detail_cigar.hybrid_length, resize_ratio))):
                # print(len(base_cigar_maps), i, ref_pointer_for_map, ref_pointer_for_draw)

                draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], None), ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], "D"), ref_pointer_for_draw, read_pointer_for_draw + gap_width, "down")

                ref_pointer_for_map = (detail_cigar.ref_start - left_extension_start) + int(i * resize_ratio)

                ref_pointer_for_draw += 1

        elif detail_cigar.op in ["I", "MappedI", "UnmappedI"]:
            # # draw duplication line and cigar map
            if detail_cigar.op == "MappedI":
                source_ref_start, source_ref_end, source_strand = detail_cigar.secondary_ref_start, detail_cigar.secondary_ref_end, detail_cigar.secondary_ref_strand
                source_ref_span = source_ref_end - source_ref_start

                # # draw cigar map for duplication
                ref_pointer_for_map = detail_cigar.ref_start - left_extension_start
                ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

                source_ref_pointer_for_map = source_ref_start - left_extension_start
                source_ref_pointer_for_draw = int(source_ref_pointer_for_map / resize_ratio)

                line_length = 0

                # # NOTE: in case that the line is out of the scope of target img
                target_img_length_on_ref_axis = np.shape(target_img)[1]
                resized_source_ref_span = int((source_ref_end - source_ref_start) / resize_ratio)
                if source_ref_pointer_for_draw < 0:
                    source_ref_pointer_for_draw = 0
                if source_ref_pointer_for_map < 0:
                    source_ref_pointer_for_map = 0
                if source_ref_pointer_for_draw + resized_source_ref_span >= target_img_length_on_ref_axis:
                    source_ref_pointer_for_draw = target_img_length_on_ref_axis - resized_source_ref_span - 1

                if source_ref_pointer_for_map + source_ref_span >= len(base_cigar_maps):
                    source_ref_pointer_for_map = (right_extension_start + right_extension_length - source_ref_span) - left_extension_start - 1

                # # set cigar's plot info
                detail_cigar.set_plot_info(source_ref_pointer_for_draw, read_pointer_for_draw, resize_ratio)

                # # draw
                if source_strand == "+":
                    for i in np.arange(0, source_ref_span, resize_ratio):
                        # # the base map if also a duplication

                        if base_cigar_maps[source_ref_pointer_for_map].freq_map["I"] != 0:
                            draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], "I", mappedI="mapped"), source_ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                        else:
                            draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[source_ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], None), source_ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                        draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], "I", mappedI="mapped"), source_ref_pointer_for_draw, read_pointer_for_draw + gap_width, "down")

                        source_ref_pointer_for_map += int(resize_ratio)
                        source_ref_pointer_for_draw += 1

                        read_pointer_for_draw += 1

                        line_length += 1
                else:
                    for i in np.arange(0, source_ref_span, resize_ratio):
                        read_pointer_for_draw += 1

                    for i in np.arange(0, source_ref_span, resize_ratio):

                        # # the base map if also a duplication
                        # if base_cigar_maps[ref_pointer_for_map].freq_map["I"] != 0 and base_cigar_maps[source_ref_pointer_for_map].freq_map["I"] != 0:
                        if base_cigar_maps[source_ref_pointer_for_map].freq_map["I"] != 0:
                            draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], "I", mappedI="inv_mapped"), source_ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                        else:
                            draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[source_ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], None), source_ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                        draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[source_ref_pointer_for_map], "I", mappedI="inv_mapped"), source_ref_pointer_for_draw, read_pointer_for_draw + gap_width, "down")

                        source_ref_pointer_for_map += int(resize_ratio)
                        source_ref_pointer_for_draw += 1
                        read_pointer_for_draw -= 1

                        line_length += 1

                    for i in np.arange(0, source_ref_span, resize_ratio):
                        read_pointer_for_draw += 1

                lines.append([source_ref_pointer_for_draw - line_length, source_ref_pointer_for_draw, read_pointer_for_draw - line_length, read_pointer_for_draw, source_strand])

            # # draw inserted seq
            else:
                ref_pointer_for_map = detail_cigar.ref_start - left_extension_start
                ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

                # # set cigar's plot info
                detail_cigar.set_plot_info(ref_pointer_for_draw, read_pointer_for_draw, resize_ratio)

                inserted_seq = detail_cigar.alt_bases

                for i in np.arange(0, len(inserted_seq), resize_ratio):
                    inserted_base = inserted_seq[int(i)]

                    # # base is also "I"
                    if base_cigar_maps[ref_pointer_for_map].freq_map["I"] != 0:
                        base_inserted_seq = base_cigar_maps[ref_pointer_for_map].map["I"][0].alt_bases

                        if i < len(base_inserted_seq):
                            base_inserted_base = base_inserted_seq[int(i)]
                            
                            draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], "I", inserted_base=inserted_base), ref_pointer_for_draw - gap_width, read_pointer_for_draw, "left")
                    
                    else:
                        draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], None), ref_pointer_for_draw - gap_width, read_pointer_for_draw, "left")

                    draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], "I", inserted_base=inserted_base), ref_pointer_for_draw + gap_width, read_pointer_for_draw, "right")

                    read_pointer_for_draw += 1

        elif detail_cigar.op == "V":

            ref_pointer_for_map = detail_cigar.ref_start - left_extension_start
            ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

            line_length = 0

            # # set cigar's plot info
            detail_cigar.set_plot_info(ref_pointer_for_draw, read_pointer_for_draw, resize_ratio)

            # # why add op len first? since we need adjust the cords when meeting reversed cigar
            for i in range(len(np.arange(0, detail_cigar.hybrid_length, resize_ratio))):
                read_pointer_for_draw += 1

            for i in range(len(np.arange(0, detail_cigar.hybrid_length, resize_ratio))):

                draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], None), ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], "V"), ref_pointer_for_draw, read_pointer_for_draw + gap_width, "down")

                ref_pointer_for_map = (detail_cigar.ref_start - left_extension_start) + int(i * resize_ratio)
                ref_pointer_for_draw += 1
                read_pointer_for_draw -= 1

                line_length += 1

            # # re-add read pointer since we subtract it in the above iteration
            for i in range(len(np.arange(0, detail_cigar.hybrid_length, resize_ratio))):
                read_pointer_for_draw += 1

            lines.append([ref_pointer_for_draw - line_length, ref_pointer_for_draw, read_pointer_for_draw - line_length, read_pointer_for_draw, "-"])

        else:
            ref_pointer_for_map = detail_cigar.ref_start - left_extension_start
            ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

            # # set cigar's plot info
            detail_cigar.set_plot_info(ref_pointer_for_draw, read_pointer_for_draw, resize_ratio)

            for i in range(len(np.arange(0, detail_cigar.hybrid_length, resize_ratio))):
                draw_single_cigar_map(target_img, get_drawable_map(base_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], None), ref_pointer_for_draw, read_pointer_for_draw - gap_width, "up")
                draw_single_cigar_map(target_img, get_drawable_map(target_cigar_maps[ref_pointer_for_map], ref_cigar_maps[ref_pointer_for_map], detail_cigar.op), ref_pointer_for_draw, read_pointer_for_draw + gap_width, "down")

                ref_pointer_for_map = (detail_cigar.ref_start - left_extension_start) + int(i * resize_ratio)
                ref_pointer_for_draw += 1
                read_pointer_for_draw += 1

    # # STEP: draw part 3, draw right extensions, including both line and maps
    ref_pointer_for_map = right_extension_start - left_extension_start
    ref_pointer_for_draw = int(ref_pointer_for_map / resize_ratio)

    line_length = 0

    for i in np.arange(0, right_extension_length, resize_ratio):
        # target_img[max(read_pointer - gap_length - bar_width, 0): read_pointer - gap_length, ref_pointer, :] = [227, 227, 227]
        # target_img[read_pointer + gap_length: read_pointer + gap_length + bar_width, ref_pointer, :] = [227, 227, 227]

        # # why we do not use slice assignment (the two lines above) for the target img: to avoid that the right extension might cover the non-match region
        # # so we use == [255, 255, 255] to determine if this pos has already been drawn
        for tmp_read_pointer in range(max(read_pointer_for_draw - gap_width - bar_width, 0), read_pointer_for_draw - gap_width):
            if tmp_read_pointer >= img_length_on_read_axis or ref_pointer_for_draw >= img_length_on_ref_axis:
                continue
            if list(target_img[tmp_read_pointer, ref_pointer_for_draw, :]) == [255, 255, 255]:
                target_img[tmp_read_pointer, ref_pointer_for_draw, :] = OP_COLOR_DICT["N"]

        for tmp_read_pointer in range(read_pointer_for_draw + gap_width, min(read_pointer_for_draw + gap_width + bar_width, img_length_on_read_axis)):
            if tmp_read_pointer >= img_length_on_read_axis or ref_pointer_for_draw >= img_length_on_ref_axis:
                continue
            if list(target_img[tmp_read_pointer, ref_pointer_for_draw, :]) == [255, 255, 255]:
                target_img[tmp_read_pointer, ref_pointer_for_draw, :] = OP_COLOR_DICT["N"]

        ref_pointer_for_draw += 1
        read_pointer_for_draw += 1

        line_length += 1

    lines.append([ref_pointer_for_draw - line_length, ref_pointer_for_draw, read_pointer_for_draw - line_length, read_pointer_for_draw, "+"])


    # # draw those lines into target img
    for line_ref_start, line_ref_end, line_read_start, line_read_end, line_strand in lines:
        draw_single_line(target_img, line_ref_start, line_ref_end, line_read_start, line_read_end, line_strand)


def calculate_img_size(left_extension_length, right_extension_length, detail_cigars):
    """
    calculate img size by cigar and extensions
    :param left_extension_length: 
    :param right_extension_length: 
    :param included_non_match_cigars: 
    :return: 
    """

    # # init
    ref_length = 0
    read_length = 0

    # # add left extension length
    ref_length += left_extension_length
    read_length += left_extension_length

    # # calculate ref and read length by cigar map
    for cigar in detail_cigars:
        if cigar.op == "D":
            ref_length += cigar.hybrid_length
        elif cigar.op in ["I", "MappedI", "UnmappedI"]:
            insert_length = cigar.hybrid_length
            read_length += insert_length
        else:
            ref_length += cigar.hybrid_length
            read_length += cigar.hybrid_length

    # # add the right extension length
    ref_length += right_extension_length
    read_length += right_extension_length

    return [read_length, ref_length]


def draw_rectangle_bounding_box(target_img, ref_pointer, read_pointer, var_ref_length, var_read_length, direct="horizontal"):
    """
    draw bounding box

    :return:
    """

    img_length_on_read_axis, img_length_on_ref_axis, _ = np.shape(target_img)

    gap_length = 2

    bar_width = int(img_length_on_read_axis / 10)

    # bounding_box_extra_extension = int(var_ref_length * 0.1)
    bounding_box_extra_extension = 10
    if direct == "horizontal":
        bounding_box_ref_start = ref_pointer - bounding_box_extra_extension
        bounding_box_read_start = read_pointer - bounding_box_extra_extension - bar_width - gap_length

        bounding_box_ref_span = 2 * bounding_box_extra_extension + var_ref_length
        bounding_box_read_span = 2 * bounding_box_extra_extension + var_read_length + 2 * bar_width + 2 * gap_length

        bounding_box_ref_end = bounding_box_ref_start + bounding_box_ref_span
        bounding_box_read_end = bounding_box_read_start + bounding_box_read_span

        cv.rectangle(target_img, (bounding_box_ref_start, bounding_box_read_start), (bounding_box_ref_end, bounding_box_read_end), [0, 0, 255], 3)

        return [(bounding_box_ref_start, bounding_box_read_start), (bounding_box_ref_end, bounding_box_read_end)]

    else:
        bounding_box_ref_start = ref_pointer - bounding_box_extra_extension - bar_width - gap_length
        bounding_box_read_start = read_pointer - bounding_box_extra_extension

        bounding_box_ref_span = 2 * bounding_box_extra_extension + var_ref_length + 2 * bar_width + 2 * gap_length
        bounding_box_read_span = 2 * bounding_box_extra_extension + var_read_length

        bounding_box_ref_end = bounding_box_ref_start + bounding_box_ref_span
        bounding_box_read_end = bounding_box_read_start + bounding_box_read_span

        cv.rectangle(target_img, (bounding_box_ref_start, bounding_box_read_start), (bounding_box_ref_end, bounding_box_read_end), [0, 0, 255], 3)

        return [(bounding_box_ref_start, bounding_box_read_start), (bounding_box_ref_end, bounding_box_read_end)]


def generate_color_plot_for_partition(partition, partial_img_path, options):
    """
    generate cigar map img for both base and target file (normal and tumor file)
    :return:
    """

    # # STEP: adjust draw parameters according to imgs size
    global line_width, gap_width, bar_width
    if options.img_size == 256:
        line_width = 1
        gap_width = 1
        bar_width = 25
    elif options.img_size == 512:
        line_width = 2
        gap_width = 2
        bar_width = 50
    elif options.img_size == 1024:
        line_width = 4
        gap_width = 4
        bar_width = 100
    else:
        pass

    detail_cigars = partition.get_detail_cigars()

    target_cigar_maps = partition.target_cigar_maps
    ref_cigar_maps = partition.ref_cigar_maps

    out_img_file_names = []
    for base_file_index in partition.base_cigar_maps.keys():
        base_cigar_maps = partition.base_cigar_maps[base_file_index]

        # # STEP: calculate extension length
        left_extension_start, right_extension_end = partition.left_extension_start, partition.right_extension_end
        left_extension_end, right_extension_start = partition.ref_start - 1, partition.ref_end + 1

        partition_ref_start, partition_ref_end = partition.ref_start, partition.ref_end

        left_extension_length = partition_ref_start - left_extension_start
        right_extension_length = right_extension_end - partition_ref_end

        # # STEP: init img, generate empty img and calculate resize ratio when given the specific output img size
        img_length_on_read_axis, img_length_on_ref_axis = calculate_img_size(left_extension_length, right_extension_length, detail_cigars)
        resize_ratio = max(img_length_on_ref_axis, img_length_on_read_axis) / options.img_size

        # print(img_length_on_ref_axis, img_length_on_read_axis, resize_ratio)
        # print(resize_ratio, img_length_on_read_axis, img_length_on_ref_axis)
        partition.set_color_plot_resize_ratio(resize_ratio)

        color_plot_img = np.ones((options.img_size, options.img_size, 3)) * 255

        # # STEP: draw lines and cigar maps
        draw_lines_and_cigar_maps(color_plot_img, left_extension_start, right_extension_start, left_extension_length, right_extension_length, detail_cigars, target_cigar_maps, base_cigar_maps, ref_cigar_maps, resize_ratio)

        out_img_file_name = "{}.{}_{}_{}_{}_{}_base{}.png".format(options.sample_name, partition.ref_chrom, partition.ref_start, partition.ref_end, partition.target_variant_type, partition.hybrid_length, base_file_index)

        cv.imwrite(os.path.join(partial_img_path, out_img_file_name), color_plot_img)

        out_img_file_names.append(out_img_file_name)

    return out_img_file_names

