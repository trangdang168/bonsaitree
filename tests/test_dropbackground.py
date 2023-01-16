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

    toy_dict2 = {
        '4': [1, None, 'd', 'e'],
        '5': [2, None, 'd', 'e'],
        'e': [2, None, 'l', 'm'],
    }
