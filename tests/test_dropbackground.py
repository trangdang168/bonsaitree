import copy
import gzip
import json
import logging
import os
from collections import defaultdict
from itertools import combinations, product

import numpy as np
import pytest
import random

from ..bonsaitree import utils
from ..bonsaitree import distributions
from ..bonsaitree import point_predictor

from ..bonsaitree.bonsai import build_pedigree
from ..bonsaitree.copytools import deepcopy
from ..bonsaitree.pedigree_object import PedigreeObject, merge_pedigrees_on_founder, extend_up, extend_down, connect_pedigrees_through_founders
from ..bonsaitree.exceptions import *
from ..bonsaitree.analytical_likelihood_tools import *
from ..bonsaitree.connect_pedigree_tools import *
from ..bonsaitree.node_dict_tools import *
from ..bonsaitree.ibd_tools import *
from ..bonsaitree.build_pedigree import *
from ..bonsaitree.distributions import AUTO_GENOME_LENGTH as GENOME_LENGTH

from .convert_ids import translate_ids_to_nums_pedigrees, translate_ids_to_nums_ibds

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")

INF = float('inf')

log = logging.getLogger(__name__)

# ================== Test utils

# goal: accurately remove the background ibd according to the current assumption
"""
INPUT FORMATS

IBD format:  genotype_id_1, genotype_id_2, chromosome, start, end, is_full_ibd, length (cM)
Pedigree Dict: the key is the id, the value is [sex, age, parent_1, parent_2]
Sex: Male = 1, Female = 0, None = None

"""


def test_one():
    up_dict1 = {
        1 : [None, None, -1, -2],
        2 : [None, None, -1, -2],
        3 : [None, None, -3, -4],
        4 : [None, None, -3, -4],
        -1 : [None, None, -5, -6],
        -3 : [None, None, -5, -7],
    }
    po1 = PedigreeObject(up_dict1)

    up_dict2 = {
        5 : [None, None]
    }
    po2 = PedigreeObject(up_dict2)

    ibd_seg_list = [
        [1, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '17', 12456159, 81041077, False, 96.75493621826172],
        [4, 5, '17', 55701470, 80890638, False, 46.832611083984375],
    ]

    new_ca1 = drop_background_ibd(
        ca1 = -5,
        ca2 = 5,
        po1 = po1,
        po2 = po2,
        ibd_seg_list = ibd_seg_list,
        alpha = 0.01,
    )

    assert(new_ca1 == -3)

def test_two():
    up_dict1 = {
        1 : [None, None, -1, -2],
        2 : [None, None, -1, -2],
        3 : [None, None, -3, -4],
        4 : [None, None, -3, -4],
        -1 : [None, None, -5, -6],
        -3 : [None, None, -5, -7],
    }
    po1 = PedigreeObject(up_dict1)

    up_dict2 = {
        5 : [None, None]
    }
    po2 = PedigreeObject(up_dict2)

    ibd_seg_list = [
        [1, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '17', 12456159, 81041077, False, 96.75493621826172],
        [4, 5, '17', 37706267, 66391837, False, 28.025470733642578],
    ]

    new_ca1 = drop_background_ibd(
        ca1 = -5,
        ca2 = 5,
        po1 = po1,
        po2 = po2,
        ibd_seg_list = ibd_seg_list,
        alpha = 0.01,
    )

    assert(new_ca1 == 3)


def test_three():
    up_dict1 = {
        1 : [None, None, -1, -2],
        2 : [None, None, -1, -2],
        3 : [None, None, -3, -4],
        4 : [None, None, -3, -4],
        -1 : [None, None, -5, -6],
        -3 : [None, None, -5, -7],
    }
    po1 = PedigreeObject(up_dict1)

    up_dict2 = {
        5 : [None, None]
    }
    po2 = PedigreeObject(up_dict2)

    ibd_seg_list = [
        [1, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '3', 37706267, 66391837, False, 28.025470733642578],
        [3, 5, '17', 55701470, 80890638, False, 46.832611083984375],
        [4, 5, '17', 55701470, 80890638, False, 46.832611083984375],
    ]

    new_ca1 = drop_background_ibd(
        ca1 = -5,
        ca2 = 5,
        po1 = po1,
        po2 = po2,
        ibd_seg_list = ibd_seg_list,
        alpha = 0.01,
    )

    assert(new_ca1 == -5)

def test_toy1():
    toy_dict1 = {
                    '1': [0, None, 'b', 'a'],
                    '2': [1, None, 'b', 'a'],
                    '3': [1, None, 'b', 'a'],
                    'a': [0, None, '0', '0'],
                    'b': [1, None, 'h', 'g'],
                    'c': [0, None, 'h', 'g'],
                    # 'd': [1, None, 'h', 'g'],
                    'g': [0, None, '0', '0'],
                    'h': [1, None, 'l', 'k'],
                    'i': [0, None, 'n', 'o'],
                    'j': [1, None, '0', '0'],
                    'k': [0, None, '0', '0'],
                    'l': [1, None, 'p', 'q'],
                }

    po1 = PedigreeObject(translate_ids_to_nums_pedigrees(toy_dict1))

    toy_dict2 = {
        '4': [1, None, 'd', 'e'],
        '5': [2, None, 'd', 'e'],
        'e': [2, None, 'l', 'm'],
    }

    po2 = PedigreeObject(translate_ids_to_nums_pedigrees(toy_dict2))

    # original segments
    ibd_seg_list_truth =[
        ['1', '2', '21', 15225121, 31491707, False, 12.383539],
        # ['1', '2', '21', 10867286, 11158632, False, .911326],
        ['1', '3', '21', 11159780, 31491707, False, 14.737657],
        ['1', '3', '21', 11159780, 31491707, False, 14.737657],
        ['1', '3', '21', 31492197, 48097120, False, 15.100822],
        # ['1', '3', '21', 10867286, 11158632, False, .911326],
        ['1', '5', '21', 10867286, 48097120, False, 30.747429],
        ['2', '3', '21', 15225121, 48097120, False, 27.481983],
        ['2', '6', '21', 10867286, 36226702, False, 27.282493],
        ['2', '6', '21', 42427818, 48097120, False, 2.299955],
        ['2', '7', '21', 10867286, 36226702, False, 27.282493],
        ['2', '7', '21', 42427818, 48097120, False, 2.299955],
        ['2', '8', '21', 10867286, 36226702, False, 27.282493],
        ['2', '8', '21', 42427818, 48097120, False, 2.299955],
        ['2', '4', '21', 42427818, 48097120, False, 2.299955],
        # ['2', '4', '21', 36228656, 42336643, False, -2.1892]
        ['2', '5', '21', 37477939, 42427353, False, 3.675311],
        ['3', '5', '21', 10867286, 48097120, False, 30.747429],
        ['6', '7', '21', 10867286, 48097120, False, 30.747429],
        ['6', '8', '21', 42186395, 48097120, False, -3.577975],
        ['6', '8', '21', 38091481, 42182870, False, 8.217810],
        ['6', '8', '21', 38091481, 42182870, False, 8.217810],
        ['6', '8', '21', 10867286, 38090081, False, 26.106673],
        ['4', '6', '21', 42342787, 48097120, False, 5.659791],
        ['7', '8', '21', 42186395, 48097120, False, -3.577975],
        ['7', '8', '21', 42186395, 48097120, False, -3.577975],
        ['7', '8', '21', 38091481, 42182870, False, 8.217810],
        ['7', '8', '21', 10867286, 38090081, False, 26.106673],
        ['7', '8', '21', 10867286, 38090081, False, 26.106673],
        ['4', '7', '21', 42342787, 48097120, False, 5.659791],
        ['4', '8', '21', 42342787, 48097120, False, 5.659791],
        ['4', '4', '21', 39834815, 41443417, False, -8.113669],
        ['4', '5', '21', 22817497, 42336643, False, 15.385181],
]

    ibd_seg_list_germline = [
        ['1', '4',	'21', 15319116, 21839933, False, 15.194	],
        ['2', '4',	'21', 15319116, 21839933, False, 15.194	],
        ['3', '4',	'21', 15319116, 21839933, False, 15.194	],
        ['4', '4',	'21', 10867286, 21839933, False, 16.760	],
        ['1', '4',	'21', 15319116, 22782845, False, 17.214	],
        ['2', '4',	'21', 15319116, 22782845, False, 17.214	],
        ['3', '4',	'21', 15319116, 22782845, False, 17.214	],
        ['4', '5',	'21', 22782861, 29985905, False, 10.909	],
        ['1', '2',	'21', 15319116, 31427808, False, 	29.053],	
        ['1', '3',	'21', 11297109, 31427808, False, 	30.088],	
        ['2', '6',	'21', 10867286, 36173761, False, 	37.559],	
        ['2', '7',	'21', 10867286, 36173761, False, 	37.559],	
        ['2', '8',	'21', 10867286, 36173761, False, 	37.559],	
        ['4', '5',	'21', 30124343, 37478248, False, 10.308	],
        ['7', '8',	'21', 10867286, 38039372, False, 	40.986],	
        ['2', '4',	'21', 39856727, 41396998, False, 3.589],	
        ['4', '4',	'21', 39856727, 41396998, False, 3.589],	
       ['4', '5',	'21', 39856727, 41396998, False, 3.589],	
       ['6', '8',	'21', 38181451, 42171704, False, 7.848],	
        ['2', '4',	'21', 36317666, 42317996, False, 12.016	],
        ['2', '5',	'21', 37480993, 42317996, False, 10.151	],
        ['4', '5',	'21', 37480993, 42317996, False, 10.151	],
        ['1', '1',	'21', 43872564, 46307297, False, 4.310],        
        ['1', '3',	'21', 43872564, 46307297, False, 4.310],
        ['1', '5',	'21', 43872564, 46307297, False, 4.310],        
        ['1', '3',	'21', 10867286, 48097120, False, 	62.728],	
        ['1', '5',	'21', 10867286, 48097120, False, 	62.728],	
        ['2', '6',	'21', 42490623, 48097120, False, 12.133	],
        ['2', '7',	'21', 42490623, 48097120, False, 12.133	],
        ['2', '8',	'21', 42490623, 48097120, False, 12.133	],
        ['2', '4',	'21', 42490623, 48097120, False, 12.133	],
        ['2', '3',	'21', 15319116, 48097120, False, 	61.162],	
        ['3', '5',	'21', 10867286, 48097120, False, 	62.728],	
        ['6', '7',	'21', 10867286, 48097120, False, 	62.728],	
        ['6', '8',	'21', 10867286, 48097120, False, 	62.728],	
        ['6', '4',	'21', 42490623, 48097120, False, 12.133	],
        ['7', '8',	'21', 10867286, 48097120, False, 	62.728],	
        ['7', '4',	'21', 42490623, 48097120, False, 12.133	],
        ['7', '8',	'21', 42172232, 48097120, False, 13.294	],
        ['8', '4',	'21', 42490623, 48097120, False, 12.133	],
    ]

    ibd_seg_list_germline = translate_ids_to_nums_ibds(ibd_seg_list_germline)

    new_ca1 = drop_background_ibd(
        ca1 = -4,
        ca2 = -5,
        po1 = po1,
        po2 = po2,
        ibd_seg_list = ibd_seg_list_germline,
        alpha = 0.01,
    )

    assert(new_ca1 == -13)