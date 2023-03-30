from bonsaitree.pedigree_object import PedigreeObject
from bonsaitree.ibd_tools import get_segment_length_list,  merge_ibd_segs
from bonsaitree.node_dict_tools import extract_down_node_dict
from bonsaitree.analytical_likelihood_tools import get_log_prob_ibd, get_ibd_pattern_log_prob, get_background_test_pval_gamma, get_var_total_length_approx, get_log_like_total_length_normal, get_background_test_pval_normal
from bonsaitree.distributions import load_distributions, AUTO_GENOME_LENGTH, MIN_PARENT_CHILD_AGE_DIFF

from bonsai_tools import *
from tests.convert_ids import *
from druid_amish import *

# class IBD:
#     def __init__(self, id1, id2, chrom, start, end, is_full, cm):
#         self.id1, self.id2, self.chrom, self.start, self.end, self.is_full, self.cm = id1, id2, chrom, start, end, is_full, cm
    
#     def __str__(self):
#         return str((self.id1, self.id2, self.chrom, self.start, self.end, self.is_full, self.cm))

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
                chrom_complete_ibd_segs_dict[chrom] = []
            chrom_ibd_segs_dict[chrom].append((start,end))
            chrom_complete_ibd_segs_dict[chrom].append([id1, id2, chrom, start, end, is_full, cm])
    return chrom_ibd_segs_dict, chrom_complete_ibd_segs_dict

def group_ibd_pairs(chrom_complete_ibd_segs_dict):
    chrom_ibd_cohorts = {} # {"<chromosome> <start> <end> <is_full_ibd>"}
    IBD_set = {}

    for chrom in chrom_complete_ibd_segs_dict:
        IBD_set[chrom] = {}
        for ibd_seg in chrom_complete_ibd_segs_dict[chrom]:
            # "<chromosome> <start> <end> <is_full_ibd>"
            skeleton = tuple(ibd_seg[2:])

            # find pair of individuals
            indv1 = ibd_seg[0]
            indv2 = ibd_seg[1]

            if skeleton not in IBD_set[chrom]:
                IBD_set[chrom][skeleton] = [[indv1, indv2]]
            else:
                IBD_set[chrom][skeleton].append([indv1, indv2])
    
    for chrom in IBD_set:
        chrom_ibd_cohorts[chrom] = {}
        for skeleton, indv_pairs in IBD_set[chrom].items():
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
            for group in groups:
                # this last argument is consistent with json file
                ibd = skeleton + (min(group),)

                # We add the individuals into the IBDs in a sorted order. Currently,
                # we are adding indvs to IBD once (here), so we sort them only once.
                chrom_ibd_cohorts[chrom][ibd] = group
    return chrom_ibd_cohorts

def get_length_and_cohorts(chrom_ibd_cohorts):
    
    chrom_ibd_lengths = {}
    chrom_ibd_members = {}

    for chrom in chrom_ibd_cohorts:
        chrom_ibd_lengths[chrom] = {}
        chrom_ibd_members[chrom] = {}
        for ibd in chrom_ibd_cohorts[chrom]:
            chrom_ibd_lengths[chrom][ibd] = ibd[-2] # get cm length
            chrom_ibd_members[chrom][ibd] = chrom_ibd_cohorts[chrom][ibd]
    
    return chrom_ibd_lengths, chrom_ibd_members


def calculate_removal_orders(chrom_complete_ibd_segs_dict, indept_gt_set1, indept_gt_set2):
    chrom_ibd_cohorts = group_ibd_pairs(chrom_complete_ibd_segs_dict)
    chrom_ibd_lengths, chrom_ibd_members = get_length_and_cohorts(chrom_ibd_cohorts)
    chrom_ibd_removal_orders = {}
    for chrom in chrom_ibd_cohorts:
        chrom_ibd_removal_orders[chrom] = sorted(chrom_ibd_cohorts[chrom], 
        key = lambda ibd: (chrom_ibd_lengths[chrom][ibd], -len(chrom_ibd_members[chrom][ibd] & indept_gt_set1), -len(chrom_ibd_members[chrom][ibd] & indept_gt_set2)))
    
    return chrom_ibd_removal_orders


def estimate_ibd_prob(node_dict, ca1, ca2, indept_gt_set1, indept_gt_set2, root_id):
    log_prob_ibd = get_log_prob_ibd(
        node_dict = node_dict,
        root_id = root_id,
        left_common_anc = ca1,
        right_common_anc = ca2,
        num_common_ancs = 1,
        left_indep_leaf_set = indept_gt_set1,
        right_indep_leaf_set = indept_gt_set2,
    )
    prob_ibd = np.exp(log_prob_ibd)
    mean = prob_ibd * AUTO_GENOME_LENGTH
    var,El,El2 =  get_var_total_length_approx(
        node_dict = node_dict,
        indep_leaf_set1 = indept_gt_set1,
        indep_leaf_set2 = indept_gt_set2,
        root_id = root_id,
        left_common_anc = ca1,
        right_common_anc = ca2,
        num_common_ancs = 1,
    )

    return mean, var

def drop_background(node_dict, ca1, indept_gt_set1, ca2, indept_gt_set2, root_id, ibd_seg_list):
    alpha = 0.01
    # retrieve all ibd segments shared between two sets, and compute total length
    chrom_seg_dict, chrom_complete_seg_dict = get_ibd_segs_between_sets(indept_gt_set1, indept_gt_set2, ibd_seg_list)
    merged_chrom_seg_dict = merge_ibd_segs(chrom_seg_dict)
    merged_seg_lengths = get_segment_length_list(merged_chrom_seg_dict)
    L_tot = sum(merged_seg_lengths)

    chrom_ibd_removal_orders = calculate_removal_orders(chrom_complete_seg_dict, indept_gt_set1, indept_gt_set2)

    mean, var = estimate_ibd_prob(node_dict, ca1, ca2, indept_gt_set1, indept_gt_set2, root_id)

    pval = get_background_test_pval_gamma(L_tot,mean,var) # Get the background IBD pvalue

    removed_ibd_index = -1
    chrom = '21' # testing on only one chromosome
    while pval < alpha and removed_ibd_index < len(chrom_ibd_removal_orders[chrom]) - 1:
        removed_ibd_index +=1
        removed_ibd_length = chrom_ibd_removal_orders[chrom][removed_ibd_index][-1]
        L_tot -= removed_ibd_length
        pval = get_background_test_pval_gamma(L_tot,mean,var)
    
    return chrom_ibd_removal_orders[chrom][removed_ibd_index:]

# gathering inputs for drop background
def get_mcra_node_dict(mcra):
    node_dicts = {}
    founders = [-1774, -2080, -91862]
    node_dict_path = "/homes/thdang/trang_amish/node_dict_"
    for founder in founders:
        with open(node_dict_path + str(-int(founder)) + ".json") as json_file:
            data = json.load(json_file)
            node_dicts[founder] = data

    for founder in founders:
        founder_node_dict = node_dicts[founder]
        if mcra in founder_node_dict:
            return extract_down_node_dict(mcra, founder_node_dict)

    return None



