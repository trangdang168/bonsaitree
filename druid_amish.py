import ast
import json
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys

from itertools import permutations, product

from bonsaitree.connect_pedigree_tools import drop_background_ibd, infer_degree_generalized_druid
from bonsaitree.pedigree_object import PedigreeObject
from bonsaitree.node_dict_tools import get_node_dict, get_node_dict_for_root, extract_down_node_dict, get_desc_deg_dict, get_min_id 
from bonsaitree.ibd_tools import check_overlap, get_segment_length_list, get_ibd_segs_between_sets, merge_ibd_segs, get_related_sets

from bonsai_tools import *
from tests.convert_ids import *

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

SEPARATOR = ";"

def druid(ca1, ca2, po, ibd_seg_list):
    # if ca1 > 0:
    #     return None, ca1, 0

    ca1_desc_set = po.rel_dict[ca1]['desc']
    ca2_desc_set = po.rel_dict[ca2]['desc']

    logging.debug("ca1 desc set " + str(ca1_desc_set))
    logging.debug("ca2 desc set " + str(ca2_desc_set))

    gt_set1 = {i for i in ca1_desc_set if i > 0}
    gt_set2 = {i for i in ca2_desc_set if i > 0}

    # ca1 or ca2 might be a leaf
    if ca1 > 0:
        gt_set1 |= {ca1}
    if ca2 > 0:
        gt_set2 |= {ca2}

    # ca1 and ca2 might not be the MRCA. Adjust them to be MRCAs
    common_anc_dict1 = po.get_common_ancestor_dict([*gt_set1], get_mrcas=True)
    common_anc_dict2 = po.get_common_ancestor_dict([*gt_set2], get_mrcas=True)

    if ca1 not in common_anc_dict1:
        ca1 = [i for i in common_anc_dict1 if i in ca1_desc_set][0]
    if ca2 not in common_anc_dict2:
        ca2 = [i for i in common_anc_dict2 if i in ca2_desc_set][0]

    logging.debug("Adjusted ca1 is " + str(ca1) + " and ca2 is " + str(ca2) + ".")

    # get the node dict below ca1
    if gt_set1:
        indep_gt_set1 = po.get_independent_inds(gt_set1) # oldest ids in gt_set1 s.t. none is descended from another

        # one or no children of ca1? No information to evaluate background IBD
        logging.debug("Retrieve descendants from " + str(ca1))
        logging.debug("indep_gt_set1: " + str(indep_gt_set1) + "\n")
        if len(indep_gt_set1) < 2:
            # SHOULD NOT ENTER HERE
            return None, ca1, 0

        node_dict1 = get_node_dict(
            ped_obj = po,
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
         # SHOULD NOT ENTER HERE
        logging.debug("URGENT: ca1 has NO descendant")
        return None, ca1, 0

    if gt_set2:
        indep_gt_set2 = po.get_independent_inds(gt_set2) # oldest ids in gt_set2 s.t. none is descended from another

        if len(indep_gt_set2) == 1:
            desc_id = [*indep_gt_set2][0]
            deg = po.rels[desc_id][ca2][0]
            node_dict2 = {ca2 : {desc_id : deg}}
        else:
            node_dict2 = get_node_dict(
                ped_obj = po,
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

        # TODO Make sure that they made it here.

        chrom_seg_dict = get_ibd_segs_between_sets(desc_set, indep_gt_set2, ibd_seg_list)

        logging.debug("chrom_seg_dict: " + str(chrom_seg_dict) + "\n")

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

    logging.debug("For individuals " + str(ca1) + " and " + str(ca2) + " the total length L_tot_max is " + str(L_tot_max))

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

    if (exempt_desc_set is None or indep_gt_set2 is None):
        # SHOULD NEVER ENTER HERE EITHER
        return None, None, None

    logging.debug("druid_deg inputs: exempt_desc_set--" + str(exempt_desc_set) \
                + "\nindep_gt_set2--" + str(indep_gt_set2) +"\n" + "node dict2: " \
                + str(node_dict2) +"\nL tot max: " + str(L_tot_max) + "\n")

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

    child_deg = po.rels[exempt_child_id][ca1][0]
    adj_deg = druid_deg - child_deg
    
    adj_deg1 = int(np.floor(adj_deg/2))
    adj_deg2 = int(np.ceil(adj_deg/2))
    node_dict = {root_id : {ca1 : adj_deg1, ca2 : adj_deg2}}

    return node_dict, root_id, L_tot

def read_ca_ibds(df):
    bonsai_segments = []
    segments = df["ibd_str"] # df[df["sources"] == ca]["ibd_str"]
    for segment in segments:
        ibd_str, cohort = segment.split("-")
        chrom, start, end, cMlength, smallest_indv = ibd_str.split()
        cohort_dict =ast.literal_eval(cohort)
        indvs = list(cohort_dict.keys())

        for i in range(len(indvs) - 1):
            for j in range(i + 1, len(indvs)):
                bonsai_seg = [int(indvs[i]), int(indvs[j]), chrom, int(start), int(end), False, float(cMlength)]
                bonsai_segments.append(bonsai_seg)

    return bonsai_segments

def process_ca(ca, married_in, genotyped):
    # here, we listed both father and mother as the common ancestor. If both 
    will_examine_ca = []
    if "&" in ca:
        f, m = ca.split("&")
        if not int(f) in married_in and not int(f) in genotyped:
            will_examine_ca.append(-int(f))
        elif not int(f) in married_in and int(f) in genotyped:
            will_examine_ca.append(int(f))
        if not int(m) in married_in and not int(m) in genotyped:
            will_examine_ca.append(-int(m))
        elif not int(m) in married_in and int(m) in genotyped:
            will_examine_ca.append(int(m))
    else:
        if not int(ca) in married_in and not int(ca) in genotyped:
            will_examine_ca.append(-int(ca))
        elif not int(ca) in married_in and int(ca) in genotyped:
            will_examine_ca.append(int(ca)) 
    return will_examine_ca
                
def validate_druid(po, common_ancestor_pairs, married_in, genotyped, outputfile, mcras_df):
    # get real pedigree

    f = open(outputfile + "_log.csv", "w")

    res = {"real" : [], "druid": [], "ca1" : [], "ca2": [], "real_root_id": [], "real_adj1": [], "real_adj2": [], "fake_adj1": [], "fake_adj2":[]}
    f.write(SEPARATOR.join(["real", "druid", "ca1", "ca2", "real_root_id", "real_adj1", "real_adj2", "fake_adj1", "fake_adj2"]))
    f.write("\n")

    logging.info("Generating all pairs of sources for bonsai")

    # No need to convert ids because genotyped individuals have the same id in bonsai and thread
    ibd_seg_list = read_ca_ibds(mcras_df)

    for indv1, indv2 in tqdm(common_ancestor_pairs):

        # convert ids because ungenotyped individuals have positive ints for ids in thread, and negative ints for ids in bonsai.
        will_examine_indv1 = process_ca(indv1, married_in, genotyped)
        will_examine_indv2 = process_ca(indv2, married_in, genotyped)

        pairs = list(product(will_examine_indv1, will_examine_indv2))
        logging.debug("for indv1-" + indv1 + "-and indv2-" + indv2 +"-we have pairs-" + str(pairs))

        # goes through the bonsai ancestor (father and mother are separated, don't understand how remarriages are handled) and compute results
        for ca1, ca2 in pairs:
            logging.debug("for ca1-" + str(ca1) + "-and indv2-" + str(ca2) +"-we have pairs-" + str(pairs))

            node_dict, fake_root_id, L_tot = druid(ca1, ca2, po, ibd_seg_list)

            if (node_dict is None):
                res["real"].append(None)
                res["druid"].append(None)
                res["ca1"].append(ca1)
                res["ca2"].append(ca2)
                res["real_root_id"].append(None)

                res["real_adj1"].append(None)
                res["real_adj2"].append(None)

                res["fake_adj1"].append(None)
                res["fake_adj2"].append(None)

                # write to file
                line = [None, None, ca1, ca2, None, None, None, None, None]
                f.write(SEPARATOR.join(str(item) for item in line))
                f.write("\n")

                break

            # get ca1, ca2 real ancestor
            real_root_dict = po.get_common_ancestor_dict([ca1, ca2], get_mrcas=True)
            fake_adj1 = node_dict[fake_root_id][ca1]
            fake_adj2 = node_dict[fake_root_id][ca2]

            for real_root_id in real_root_dict:

                # po.rels[real_root_id][ca1] returns a tuple of (up_deg, down_deg, num_common_ancestors) Down is number of meioses from root to individual.

                adj1 = po.rels[real_root_id][ca1][1]
                adj2 = po.rels[real_root_id][ca2][1]
                real =  {real_root_id : {ca1 : adj1, ca2 : adj2}}

                res["real"].append(real)
                res["druid"].append(node_dict)
                res["ca1"].append(ca1)
                res["ca2"].append(ca2)

                # remarriage case, two parents might have different results.
                res["real_root_id"].append(real_root_id)

                res["real_adj1"].append(adj1)
                res["real_adj2"].append(adj2)

                res["fake_adj1"].append(fake_adj1)
                res["fake_adj2"].append(fake_adj2)
                line = [real, node_dict, ca1, ca2, real_root_id, adj1, adj2, fake_adj1, fake_adj2]
                f.write(SEPARATOR.join(str(item) for item in line))
                f.write("\n")

    logging.info("Finish collecting data")
    f.close()

    df = pd.DataFrame(res)
    df.to_csv(outputfile + "_final.csv")
    return df
    

def amish(file_prefix):
    logging.info("Results saved with prefix: " + file_prefix)

    logging.info("Retrieving AMISh married-in and genotyped individuals")
    married_in = "/homes/thdang/trang_amish/ped-aware/implementation/amish_married_in.csv"
    with open(married_in, "r") as m:
        married_in = ast.literal_eval(m.read())

    genotyped = "/homes/thdang/trang_amish/ped-aware/implementation/genotyped.txt"

    with open(genotyped, "r") as g:
        genotyped = ast.literal_eval(g.read())

    logging.info("Gathering the best sources for each ibd cohort")
    mcra_file = "/homes/thdang/trang_amish/ped-aware/implementation/amish3_phasedibd_final_sources.csv"
    mcras_df = pd.read_csv(mcra_file)
    sources = set(mcras_df["sources"].unique())
    logging.debug("DEBUG SOURCES")
    indvs = list(permutations(sources, 2))[:10]

    logging.info("There are " + str(len(indvs)) + " pairs of mcras to go through.")

    logging.info("Initializing AMISH pedigree object")
    state_dict = {}
    with open("/homes/thdang/trang_amish/extended_state_dict.json", "r") as json_file:
        state_dict = json.loads(json_file.read())
    po = PedigreeObject.init_from_state(state_dict)

    validate_druid(po, indvs, married_in, genotyped, file_prefix, mcras_df)

def process_toy_ancestor_node(ancestor_node_str):
    ancestor_node_str = ancestor_node_str.split("&")
    res = []
    for item in ancestor_node_str:
        res.append(TOY_DICT_TO_NUMBER[item])
    return "&".join(str(-item) for item in res)

def test_toy_ancestor_node():
    logging.info("test toy ancestor node")
    res = process_toy_ancestor_node("l&k")
    assert(res == "23&24")

    res = process_toy_ancestor_node("n&o")
    assert(res == "26&27")

    logging.info("all toy test ancestor node passed!")

def process_toy_ibds(ibd_str):
    ibd, cohort = ibd_str.split("-")
    cohort_dict = ast.literal_eval(cohort)
    new_cohort = {}
    for indv in cohort_dict:
        # only genotyped individuals have ibds, and all genotyped individuals have positive ids
        bonsai_id = TOY_DICT_TO_NUMBER[indv]
        new_cohort[bonsai_id] = cohort_dict[indv]
    return ibd + "-" + str(new_cohort)

def toy(file_prefix):
    logging.info("Test validate druid on toy")

    # include all individuals but p, q
    toy_dict1 = {
                    '1': [0, None, 'b', 'a'],
                    '2': [1, None, 'b', 'a'],
                    '3': [1, None, 'b', 'a'],
                    'a': [0, None, '0', '0'],
                    'b': [1, None, 'h', 'g'],
                    'c': [0, None, 'h', 'g'], #so that it doesnt get confused
                    '4': [1, None, 'd', 'e'],
                    '5': [2, None, 'd', 'e'],
                    'e': [2, None, 'l', 'm'],
                    'd': [1, None, 'h', 'g'],
                    'g': [0, None, '0', '0'],
                    'h': [1, None, 'l', 'k'],
                    'j': [1, None, '0', '0'],
                    'k': [0, None, '0', '0'],
                    'l': [1, None, 'p', 'q'],
                    'n': [1, None, 'p', 'q'],
                    'o': [0, None, '0', '0'], 
                    'i': [0, None, 'n', 'o'], 
                    'f': [1, None, 'j', 'i'], 
                    '6': [0, None, 'f', 'c'], 
                    '7': [0, None, 'f', 'c'], 
                    '8': [0, None, 'f', 'c'], 
                }

    po = PedigreeObject(translate_ids_to_nums_pedigrees(toy_dict1))
    married_in =  {14, 10, 24, 12, 11, 27} # ["a", "g", "k", "j", "i", "o"]
    genotyped =  {1, 2, 3, 4, 5, 6, 7, 8}
    mcra_file = "/homes/thdang/trang_amish/ped-aware/implementation/toy3_phasedibd_final_sources.csv"
    mcra_df = pd.read_csv(mcra_file)
    # parse mcra df into bonsai ids
    mcra_df["sources"] = mcra_df["sources"].map(process_toy_ancestor_node)
    mcra_df["ibd_str"] = mcra_df["ibd_str"].map(process_toy_ibds)
    #  [a&b, c&f] -- our code picked h&g (1, 2, 4) for smallest cohort, instead of d&e "b&a", "f&c", "h&g", 
    sources = pd.Series(["b&a", "f&c", "d&e", "h&g"]).map(process_toy_ancestor_node)
    indvs = list(permutations(sources, 2))[:10] 

    logging.info("individuals examined: " + str(indvs))

    df = validate_druid(po, indvs, married_in, genotyped, file_prefix, mcra_df)
    df["ca1_thread_id"] = df["ca1"].map(TOY_DICT_TO_STRING)
    df["ca2_thread_id"] = df["ca2"].map(TOY_DICT_TO_STRING)
    df["real_root_thread_id"] = df["real_root_id"].map(TOY_DICT_TO_STRING)
    df.to_csv("toy_test.csv")

    logging.info("Finish validating druid on toy")

def amish_ancestor_descendants():

    logging.info("Retrieving AMISh married-in and genotyped individuals")
    married_in = "/homes/thdang/trang_amish/ped-aware/implementation/amish_married_in.csv"
    with open(married_in, "r") as m:
        married_in = ast.literal_eval(m.read())

    genotyped = "/homes/thdang/trang_amish/ped-aware/implementation/genotyped.txt"

    with open(genotyped, "r") as g:
        genotyped = ast.literal_eval(g.read())

    logging.info("Gathering the best sources for each ibd cohort")
    mcra_file = "/homes/thdang/trang_amish/ped-aware/implementation/amish3_phasedibd_final_sources.csv"
    mcras_df = pd.read_csv(mcra_file)
    sources = set(mcras_df["sources"].unique())
    logging.debug("DEBUG SOURCES")

    logging.info("Initializing AMISH pedigree object")
    state_dict = {}
    with open("/homes/thdang/trang_amish/extended_state_dict.json", "r") as json_file:
        state_dict = json.loads(json_file.read())

    logging.info("Process most common recent ancestors")
    # go through each ca to gather their descendants
    processed_sources = [] # processed sources are nodes such as "5000&5001" splitted into two nodes with negative integer ids -5000 and -5001 because they are genotyped. Also, skip married in.
    for source in sources:
        processed_sources.extend(process_ca(source, married_in, genotyped))

    po = PedigreeObject.init_from_state(state_dict)
    res = {"ca": [], "desc": []}
    for ca in po.rel_dict.keys():
        res["ca"].append(ca)
        desc_set = po.rel_dict[ca]['desc']
        gt_set = {i for i in desc_set if i > 0}
        indep_gt_set = po.get_independent_inds(gt_set)
        res["desc"].append(indep_gt_set)
    
    pd.DataFrame(res).to_csv("/homes/thdang/trang_amish/ancestor_descendants.csv")

def main():
    # test_toy_ancestor_node()
    file_prefix = sys.argv[1]
    # toy(file_prefix + "_toy")
    amish(file_prefix)
    # amish_ancestor_descendants()

if __name__ == "__main__":
    main()