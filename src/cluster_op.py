from scipy.cluster.hierarchy import linkage, fcluster
# from sklearn.cluster import DBSCAN


def iterative_cluster(partition, hyper_cigar, options, ignore_op=False, ignore_length=False):
    """
    perform iterative cluster, given a partition and a target hyper_cigar, determine if the hyper_cigar belongs to this partition
    """

    # # different op
    # if "D" in hyper_cigar.op and "D" not in partition.included_hyper_op:
    #     return False, -1
    # if "D" in partition.included_hyper_op and "D" not in hyper_cigar.op:
    #     return False, -1

    if ignore_op is False and partition.included_hyper_op != hyper_cigar.op:
        return False, -1
    if ignore_op is True and ((partition.included_hyper_op == "I" and hyper_cigar.op == "D") or (partition.included_hyper_op == "D" and hyper_cigar.op == "I")):
        return False, -1
    # # different chrom
    if partition.ref_chrom != hyper_cigar.ref_chrom:
        return False, -1
    # # different size
    if ignore_length is False and partition.hybrid_length > hyper_cigar.hybrid_length and (hyper_cigar.hybrid_length / partition.hybrid_length) < options.size_sim_ratio:
        return False, -1
    if ignore_length is False and partition.hybrid_length < hyper_cigar.hybrid_length and (partition.hybrid_length / hyper_cigar.hybrid_length) < options.size_sim_ratio:
        return False, -1
    # # different pos
    if abs(hyper_cigar.ref_start - partition.ref_start) > options.dist_diff_length:
        return False, -1

    # # STEP: calculate size sim ratio
    if partition.hybrid_length > hyper_cigar.hybrid_length:
        size_sim = hyper_cigar.hybrid_length / partition.hybrid_length
    else:
        size_sim = partition.hybrid_length / hyper_cigar.hybrid_length

    # # STEP: cluster single bkp events (such as B)
    if not options.skip_bnd and hyper_cigar.op == "B":
        if partition.included_hyper_cigars[0].op != "B":
            return False, -1
        if hyper_cigar.detail_cigars[0].bnd_ref_start - partition.included_hyper_cigars[0].detail_cigars[0].bnd_ref_start <= options.dist_diff_length and hyper_cigar.detail_cigars[0].bnd_ref_chrom == partition.included_hyper_cigars[0].detail_cigars[0].bnd_ref_chrom:
            return True, 1
    # # STEP: cluster multi bkps events
    else:
        return True, size_sim

    return False, -1


def hierarchical_cluster(data, cluster_max_distance=0.3):
    """
    perform hierarchical cluster
    """

    Z = linkage(data, method="average", metric=span_position_distance)
    cluster_indices = list(fcluster(Z, cluster_max_distance, criterion='distance'))

    return cluster_indices


def dbscan_cluster(data):
    pass
    # db = DBSCAN(eps=0.3, min_samples=1, metric=span_position_distance).fit(data)
    # cluster_indices = db.labels_
    #
    # return cluster_indices


def length_distance(signature1, signature2):

    return abs(signature1[0] - signature2[0])


def span_position_distance(signature1, signature2):
    distance_normalizer = signature1[2]
    span1 = signature1[1] - signature1[0]
    span2 = signature2[1] - signature2[0]
    center1 = (signature1[0] + signature1[1]) // 2
    center2 = (signature2[0] + signature2[1]) // 2
    position_distance = min(abs(signature1[0] - signature2[0]), abs(signature1[1] - signature2[1]), abs(center1 - center2)) / distance_normalizer
    span_distance = abs(span1 - span2) / max(span1, span2)

    return position_distance + span_distance
