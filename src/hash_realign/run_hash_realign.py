#!/usr/bin/env python3

from src.hash_realign.classes import Sequence
from src.hash_realign.classes import Segment
from src.hash_realign.hash_aligner import HashAligner
from src.cigar_op import fetch_ref_seq, DetailCigar, HyperCigar


def select_longest(segments):
    longest_true = []
    longest_false = []
    for seg in segments:
        if seg.forward == True:
            if len(longest_true) == 0:
                longest_true = [seg]
            else:
                if abs(seg.read_end - seg.read_start) > abs(longest_true[0].read_end - longest_true[0].read_start):
                    longest_true = [seg]
                elif abs(seg.read_end - seg.read_start) == abs(longest_true[0].read_end - longest_true[0].read_start):
                    longest_true.append(seg)

        else:
            if len(longest_false) == 0:
                longest_false = [seg]
            else:
                if abs(seg.read_end - seg.read_start) > abs(longest_false[0].read_end - longest_false[0].read_start):
                    longest_false = [seg]
                elif abs(seg.read_end - seg.read_start) == abs(longest_false[0].read_end - longest_false[0].read_start):
                    longest_false.append(seg)
    after_select = []
    after_select.extend(longest_true)
    after_select.extend(longest_false)

    return after_select


def hash_realign_unmapped(partition, ref_file, k=17, min_accept=30):
    """
    Hashtable for realign unmapped seq
    :param ref: ref seq
    :param seq: unmapped seq
    :param k: k-mer size
    :param min_accept: min accept length for realigning
    :return: 
    """

    ref_seq = fetch_ref_seq(partition.ref_chrom, partition.left_extension_start, partition.right_extension_end, ref_file)

    new_detail_cigars = []
    modified_flag = False

    for detail_cigar in partition.included_hyper_cigars[0].detail_cigars:
        if detail_cigar.op != "UnmappedI":
            new_detail_cigars.append(detail_cigar)

        else:
            alt_seq = detail_cigar.alt_bases

            # # generate seq object
            ref = Sequence(ref_seq)
            read = Sequence(alt_seq)

            # # align ref2ref
            repeat_thresh = 2
            aligner_ref2ref = HashAligner(k, min_accept, 0, repeat_thresh)
            aligner_ref2ref.run(ref, ref)
            diff_segs = aligner_ref2ref.getSelfDiffSegs()

            y_hashvalues = aligner_ref2ref.hashvalues
            avoid_mers = aligner_ref2ref.avoid_kmers

            # # align ref2read
            aligner_ref2read = HashAligner(k, min_accept, 0, repeat_thresh)
            aligner_ref2read.run(read, ref, diff_segs, y_hashvalues, avoid_mers)
            # aligner_ref2read.run(read, ref)

            segments_merge = aligner_ref2read.getMergeSegments()

            # # adjust reversed segments' start and end cord
            for seg in segments_merge:
                if seg.forward is False:
                    tmp = seg.read_end
                    seg.read_end = seg.read_start
                    seg.read_start = tmp

            # # convert aligned segments to detail cigars
            if len(segments_merge) == 0:
                continue

            modified_flag = True

            segments_merge = sorted(segments_merge, key=lambda x:x.read_start)

            # # check the unmapped seq before the first segment
            if segments_merge[0].read_start - 0 > min_accept:
                new_detail_cigar = DetailCigar("UnmappedI", segments_merge[0].read_start - 0, detail_cigar.ref_chrom, detail_cigar.ref_start, detail_cigar.ref_end, detail_cigar.strand)
                new_detail_cigar.set_alt_bases(detail_cigar.alt_bases[0: segments_merge[0].read_start])

                new_detail_cigars.append(new_detail_cigar)

            for i in range(len(segments_merge)):
                cur_segment = segments_merge[i]

                # # generate detail cigar object
                new_detail_cigar = DetailCigar("MappedI", cur_segment.read_end - cur_segment.read_start + 1, detail_cigar.ref_chrom, detail_cigar.ref_start, detail_cigar.ref_end, detail_cigar.strand)
                new_detail_cigar.set_alt_bases(detail_cigar.alt_bases[cur_segment.read_start: cur_segment.read_end + 1])
                new_detail_cigar.set_secondary_mapping(detail_cigar.ref_chrom, cur_segment.ref_start + partition.left_extension_start, cur_segment.ref_end + partition.left_extension_start, "+" if cur_segment.forward else "-")

                new_detail_cigars.append(new_detail_cigar)

                # # not the last align, then we need to check the unmapped seq between cur and next segment
                if i != len(segments_merge) - 1:
                    next_segment = segments_merge[i + 1]

                    if next_segment.read_start - cur_segment.read_end > min_accept:
                        new_detail_cigar = DetailCigar("UnmappedI", next_segment.read_start - cur_segment.read_end, detail_cigar.ref_chrom, detail_cigar.ref_start, detail_cigar.ref_end, detail_cigar.strand)
                        new_detail_cigar.set_alt_bases(detail_cigar.alt_bases[cur_segment.read_end + 1: next_segment.read_start])

                        new_detail_cigars.append(new_detail_cigar)

            for seg in segments_merge:

                print(seg.to_string())
                # print(seg.read_start, seg.read_end, seg.ref_start + partition.left_extension_start, seg.ref_end + partition.left_extension_start)

            # select most long for both direction
            # if len(segments_merge) >= 2:
            #     segments_merge = select_longest(segments_merge)

    if modified_flag is True:
        partition.update_included_hyper_cigars([HyperCigar(new_detail_cigars, partition.supp_reads)])



