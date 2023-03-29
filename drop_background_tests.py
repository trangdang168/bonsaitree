from drop_background import *

# test set up 
def toy_set_up():
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

    return po, married_in, genotyped, mcra_df, sources

########### TEST SELF-WRITTEN MERGING IBDS ####################
def test_toy_get_ibd_between_sets():
    toy_po, married_in, genotyped, mcra_df, sources = toy_set_up()
    ibd_seg_list = read_ca_ibds(mcra_df)

    indept_gt_set1 = {2} 
    indep_gt_set2 = {4, 5, 6}
    chrom_seg_dict, chrom_seg_full_dict = get_ibd_segs_between_sets(indept_gt_set1, indep_gt_set2, ibd_seg_list)

    expected_chrom_seg_dict = {'21': [(10867286, 36249317),
     (42583395, 48097120), (36210416, 48097120), (37421765, 
     42425427), (15271418, 22681881), (15279616, 21882038)]}

    expected_chrom_seg_full_dict = {'21': [
        [2, 6, '21', 10867286, 36249317, False, 25.382030999999998], 
        [2, 6, '21', 42583395, 48097120, False, 5.513724999999994], 
        [4, 2, '21', 36210416, 48097120, False, 11.886703999999995], 
        [2, 5, '21', 37421765, 42425427, False, 5.0036619999999985], 
        [4, 2, '21', 15271418, 22681881, False, 7.410463], 
        [4, 2, '21', 15279616, 21882038, False, 6.602422000000001]]}
    
    assert(chrom_seg_dict == expected_chrom_seg_dict)
    assert(chrom_seg_full_dict == expected_chrom_seg_full_dict)
    logging.info("test_get_ibd_between_sets passed!")

def test_toy_group_ibd_pairs():
    chrom_complete_ibd_segs_dict = \
   {'21':
    [[2, 6, '21', 10867286, 36249317, False, 25.382030999999998],
    [2, 6, '21', 42583395, 48097120, False, 5.513724999999994],
    [4, 2, '21', 10867286, 36249317, False, 25.382030999999998],
    [2, 5, '21', 10867286, 36249317, False, 25.382030999999998], 
    [4, 2, '21', 15271418, 22681881, False, 7.410463], 
    [4, 2, '21', 15279616, 21882038, False, 6.602422000000001]
    ]}
    grouped = group_ibd_pairs(chrom_complete_ibd_segs_dict)
    assert(grouped == {'21': 
    {'21 10867286 36249317 False 25.382030999999998 2': {2, 4, 5, 6},
     '21 42583395 48097120 False 5.513724999999994 2': {2, 6},
     '21 15271418 22681881 False 7.410463 2': {2, 4}, 
     '21 15279616 21882038 False 6.602422000000001 2': {2, 4}}})

    logging.info("Test toy_group_ibd_pairs passed")

# test removal orders methods
def test_get_removal_orders():
    chrom_complete_ibd_segs_dict = {'21': [
        [2, 6, '21', 10867286, 36249317, False, 25.382030999999998], 
        [2, 6, '21', 42583395, 48097120, False, 5.513724999999994], 
        [4, 2, '21', 36210416, 48097120, False, 11.886703999999995], 
        [2, 5, '21', 37421765, 42425427, False, 5.0036619999999985], 
        [4, 2, '21', 15271418, 22681881, False, 7.410463], 
        [4, 2, '21', 15279616, 21882038, False, 6.602422000000001]
        ]}
    indept_gt_set1 = {2}
    indept_gt_set2 = {4, 5, 6}

    ibd_removal_orders = calculate_removal_orders(chrom_complete_ibd_segs_dict, indept_gt_set1, indept_gt_set2)
    expected_ibd_remove_orders = [
        '21 37421765 42425427 False 5.0036619999999985 2', 
        '21 42583395 48097120 False 5.513724999999994 2', 
        '21 15279616 21882038 False 6.602422000000001 2', 
        '21 15271418 22681881 False 7.410463 2', 
        '21 36210416 48097120 False 11.886703999999995 2', 
        '21 10867286 36249317 False 25.382030999999998 2']

    assert(ibd_removal_orders == expected_ibd_remove_orders)

    logging.info("Successfully get_remove_orders on different lengths (bp).")

    indept_gt_set1 = {2}
    indept_gt_set2 = {4, 5, 6}

    chrom_complete_ibd_segs_dict = {'21':
    [[2, 6, '21', 12583395, 18097120, False, 5.513724999999994],
    [2, 6, '21', 42583395, 48097120, False, 5.513724999999994],
    [4, 2, '21', 12583395, 18097120, False, 5.513724999999994],
    [2, 5, '21', 12583395, 18097120, False, 5.513724999999994], 
    [4, 2, '21', 15271418, 22681881, False, 7.410463], 
    [4, 2, '21', 15279616, 21882038, False, 6.602422000000001]
    ]}
    ibd_removal_orders = calculate_removal_orders(chrom_complete_ibd_segs_dict, indept_gt_set1, indept_gt_set2)
    
    expected_ibd_remove_orders = [
        '21 12583395 18097120 False 5.513724999999994 2', 
        '21 42583395 48097120 False 5.513724999999994 2', 
        '21 15279616 21882038 False 6.602422000000001 2', 
        '21 15271418 22681881 False 7.410463 2']

    assert(expected_ibd_remove_orders == ibd_removal_orders)
    logging.info("Successfully get_remove_orders on different cohort sizes (bp).")

def main():
    test_toy_group_ibd_pairs()
    # test_get_removal_orders()
    test_toy_get_ibd_between_sets()

if __name__ == "__main__":
    main()