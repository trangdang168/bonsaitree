from bonsaitree.pedigree_object import PedigreeObject
from bonsaitree.ibd_tools import get_segment_length_list,  merge_ibd_segs
from bonsaitree.analytical_likelihood_tools import get_log_prob_ibd, get_ibd_pattern_log_prob, get_background_test_pval_gamma, get_var_total_length_approx, get_log_like_total_length_normal, get_background_test_pval_normal

from bonsai_tools import *
from tests.convert_ids import *
from druid_amish import *


def get_ibd_segs_between_sets(
    leaves1,
    leaves2,
    ibd_seg_list
):
    """
    Get IBD segments between leaves1 and leaves2 and store them in a dict with chromosomes as keys
    Args:
        leaves1: First set of leaves
        leaves2: Second set of leaves
        ibd_seg_list: list of the form [[id1, id2, chromosome, start, end, is_full_ibd, seg_cm]]

    Returns 
        chrom_ibd_segs_dict: { chrom : [[start,end],[start,end],...] }
        chrom_complete_ibd_segs_dict: { chrom : [[id1, id2, chromosome, start, end, is_full_ibd, seg_cm],[id1, id2, chromosome, start, end, is_full_ibd, seg_cm],...] }
    """
    chrom_ibd_segs_dict : Dict[str, List[Tuple[Any,Any]]] = dict() # Dict of the form { chrom : [[start,end],[start,end],...] }
    chrom_complete_ibd_segs_dict = {}
    for id1, id2, chrom, start, end, is_full, cm in ibd_seg_list:
        if (id1 in leaves1 and id2 in leaves2) or (id1 in leaves2 and id2 in leaves1):
            if chrom not in chrom_ibd_segs_dict:
                chrom_ibd_segs_dict[chrom] = []
                chrom_complete_ibd_segs_dict
            chrom_ibd_segs_dict[chrom].append((start,end))
            chrom_complete_ibd_segs_dict[chrom].append([id1, id2, chrom, start, end, is_full, cm])
    return chrom_ibd_segs_dict, chrom_complete_ibd_segs_dict

def group_ibd_pairs(chrom_complete_ibd_segs_dict):
    IBD_cohorts = {} # {"<chromosome> <start> <end> <is_full_ibd>"}
    IBD_set = {}

    for chrom in chrom_complete_ibd_segs_dict:
        for ibd_seg in chrom_complete_ibd_segs_dict:
            # "<chromosome> <start> <end> <is_full_ibd>"
            skeleton = tuple(ibd_seg[2:])

            # find pair of individuals
            indv1 = ibd_seg[0]
            indv2 = ibd_seg[1]

            if skeleton not in IBD_set:
                IBD_set[chrom][skeleton] = [[indv1, indv2]]
            else:
                IBD_set[chrom][skeleton].append([indv1, indv2])
    
    for chrom in IBD_set:
        for skeleton, individual_pairs in IBD_set.items():
            groups = []
            # first go through all pairs
            for pair in indv_pairs:
                joined = []
                joined_indicies = []
                for i,group in enumerate(groups):
                    if pair[0] in group or pair[1] in group:
                        group.add(pair[0])
                        group.add(pair[1])
                        joined.append(group)
                        joined_indicies.append(i)

                # if we have joined 2 or more groups with this pair
                if len(joined) > 1:
                    new_group = set().union(*joined)
                    groups.append(new_group)
                    joined_indicies = reversed(joined_indicies)

                    # remove old groups
                    for i in joined_indicies:
                        groups.pop(i)

                # if no groups joined, create a new one from this pair
                elif len(joined) == 0:
                    groups.append({pair[0], pair[1]})
            # create an IBD instance for each group
            cohort = set()
            for group in groups:
                # this last argument is consistent with json file
                ibd = skeleton + (min(group))

                # We add the individuals into the IBDs in a sorted order. Currently,
                # we are adding indvs to IBD once (here), so we sort them only once.

                cohort.update(group)
            IBD_cohorts[chrom][ibd] = cohort
    return IBD_cohorts

def calculate_removal_orders():
    pass

def drop_background():
    pass