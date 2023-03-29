"""
Most Recent Common Ancestor

This code finds the most recent common ancestors of the *genotyped* 
descedants from the genotypes and the pedigree. 
The common ancestors are then used to check if their 
relationships are correctly computed by bonsai's generalized 
druid.
"""

from simulation_files import *

from bonsaitree.pedigree_object import PedigreeObject
from bonsaitree.node_dict_tools import get_node_dict, get_node_dict_for_root
from tests.convert_ids import *
from bonsai_tools import *

import json

def has_genotyped_descendants(po, ca):
    pass

# three founders
# p: 91862, m: 78497
# p: 1744, m: 17442
# p: 2080, m: 2081

def main():

    pedigree = amish_to_bonsai_ped()
    state_dict = {"one": "two"}
    ped_obj = PedigreeObject(pedigree)
    
    state_dict = ped_obj.get_state_dict()
    with open('extended_state_dict.json', 'w') as f:
        json.dump(state_dict, f)

    # indvs = set(pedigree.keys())
    # res = get_node_dict_for_root(-91862, ped_obj)
    # with open('node_dict_91862.txt', 'w') as f:
    #     f.write("%s\n" % str(res))

    # res1 = get_node_dict_for_root(-1744, ped_obj)
    # with open('node_dict_1744.txt', 'w') as f:
    #     f.write("%s\n" % str(res1))

    # res2 = get_node_dict_for_root(-2080, ped_obj)
    # with open('node_dict_2080.txt', 'w') as f:
    #     f.write("%s\n" % str(res2))

main()
