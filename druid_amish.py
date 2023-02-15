from bonsaitree.connect_pedigree_tools import drop_background_ibd
from bonsaitree.pedigree_object import PedigreeObject
from bonsaitree.node_dict_tools import get_node_dict, get_node_dict_for_root
from bonsai-tools import *

def druid(ca1, ca2, po, ibd_seg_list):
       if ca1 > 0:
        return ca1

    ca1_desc_set = po1.rel_dict[ca1]['desc']
    ca2_desc_set = po2.rel_dict[ca2]['desc']

    gt_set1 = {i for i in ca1_desc_set if i > 0}
    gt_set2 = {i for i in ca2_desc_set if i > 0}

    # ca1 or ca2 might be a leaf
    if ca1 > 0:
        gt_set1 |= {ca1}
    if ca2 > 0:
        gt_set2 |= {ca2}

    # ca1 and ca2 might not be the MRCA. Adjust them to be MRCAs
    common_anc_dict1 = po1.get_common_ancestor_dict([*gt_set1], get_mrcas=True)
    common_anc_dict2 = po2.get_common_ancestor_dict([*gt_set2], get_mrcas=True)

    if ca1 not in common_anc_dict1:
        ca1 = [i for i in common_anc_dict1 if i in ca1_desc_set][0]
    if ca2 not in common_anc_dict2:
        ca2 = [i for i in common_anc_dict2 if i in ca2_desc_set][0]

    # get the node dict below ca1
    if gt_set1:
        indep_gt_set1 = po1.get_independent_inds(gt_set1) # oldest ids in gt_set1 s.t. none is descended from another

        # one or no children of ca1? No information to evaluate background IBD
        if len(indep_gt_set1) < 2:
            return ca1

        node_dict1 = get_node_dict(
            ped_obj = po1,
            rels = indep_gt_set1,
        )
        node_dict1 = extract_down_node_dict(
            node = ca1,
            node_dict = node_dict1,
        )
    else:
        raise Exception("ca1 has no descendants.")


    # one or no children of ca1? No information to evaluate background IBD
    if len(node_dict1[ca1]) < 2:
        return ca1

    if gt_set2:
        indep_gt_set2 = po2.get_independent_inds(gt_set2) # oldest ids in gt_set2 s.t. none is descended from another

        if len(indep_gt_set2) == 1:
            desc_id = [*indep_gt_set2][0]
            deg = po2.rels[desc_id][ca2][0]
            node_dict2 = {ca2 : {desc_id : deg}}
        else:
            node_dict2 = get_node_dict(
                ped_obj = po2,
                rels = indep_gt_set2,
            )
            node_dict2 = extract_down_node_dict(
                node = ca2,
                node_dict = node_dict2,
            )
    elif ca2 > 0:
        indep_gt_set2 = {ca2}
        node_dict2 = {ca2 : {ca2 : 0}}
    else:
        raise Exception("ca2 has no descendants.")

    # find exempt child id
    desc_set_dict = {}
    L_tot_max = 0
    exempt_child_id = None
    exempt_desc_set = None
    exempt_node_dict = None
    for child_id in node_dict1[ca1]:
        desc_deg_dict = get_desc_deg_dict(
            anc_id = child_id,
            node_dict = node_dict1,
        )
        desc_set = {*desc_deg_dict}
        desc_set &= indep_gt_set1

        if child_id > 0:
            desc_set |= {child_id}

        desc_set_dict[child_id] = desc_set

        chrom_seg_dict = get_ibd_segs_between_sets(desc_set, indep_gt_set2, ibd_seg_list)
        merged_chrom_seg_dict = merge_ibd_segs(chrom_seg_dict)
        merged_seg_lengths = get_segment_length_list(merged_chrom_seg_dict)
        L_tot = sum(merged_seg_lengths)

        if L_tot > L_tot_max:
            L_tot_max = L_tot
            exempt_child_id = child_id
            exempt_desc_set = desc_set
            exempt_node_dict = extract_down_node_dict(
                node = exempt_child_id,
                node_dict = node_dict1,
            )

    # create a new node dict by removing exempt_child_id and
    # their descendants from node_dict1
    leave_one_out_node_dict1 = {}
    for child_id in node_dict1[ca1]:
        if child_id != exempt_child_id:
            child_node_dict = extract_down_node_dict(
                node = child_id,
                node_dict = node_dict1,
            )
            leave_one_out_node_dict1.update(child_node_dict)
    leave_one_out_node_dict1[ca1] = {i : d for i,d in node_dict1[ca1].items() if i != exempt_child_id}

    # get indep descs of ca1 who are not descs of exempt_child_id
    leave_one_out_gt_set1 = set.union(*[v for i,v in desc_set_dict.items() if i != exempt_child_id])

    # get L_tot shared between leave_one_out_gt_set1 and leaf_set2
    chrom_seg_dict = get_ibd_segs_between_sets(leave_one_out_gt_set1, indep_gt_set2, ibd_seg_list)
    merged_chrom_seg_dict = merge_ibd_segs(chrom_seg_dict)
    merged_seg_lengths = get_segment_length_list(merged_chrom_seg_dict)
    L_tot = sum(merged_seg_lengths)

    druid_deg = infer_degree_generalized_druid(
        leaf_set1 = exempt_desc_set,
        leaf_set2 = indep_gt_set2,
        node_dict1 = exempt_node_dict,
        node_dict2 = node_dict2,
        L_merged_tot = L_tot_max,
    )

    # create a node dict that approximates the degree between ca1 and ca2
    min_id1 = get_min_id(leave_one_out_node_dict1)
    min_id2 = get_min_id(node_dict2)
    min_id = min(min_id1, min_id2) - 1
    root_id = min(-1, min_id)  # id lower than any value in either node_dict
    
    adj_deg1 = int(np.floor(adj_deg/2))
    adj_deg2 = int(np.ceil(adj_deg/2))
    node_dict = {root_id : {ca1 : adj_deg1, ca2 : adj_deg2}}

    return node_dict, root_id, L_tot

def validate_druid():
    # get real pedigree
    pedigree = amish_to_bonsai_ped()
    po = PedigreeObject(pedigree)

    indvs = list(pedigree.keys())
    for i in range(len(indvs)):
        for j in range(i+1, len(indvs)):
            ca1 = indvs[i]
            ca2 = indvs[j]

            node_dict, root_id, L_tot = druid(ca1, ca2, po, ibd_seg_list)

            true_root = po.get_common_ancestor_dict([*gt_set1], get_mrcas=True)
            print(len(true_root))

            # check if true root matches the distance in node_dict