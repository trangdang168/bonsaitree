"""
Translate from toy string ids to numbers. Follow the bonsai format (positive numbers
for genotyped individuals, negative for ungenotyped)
"""
TOY_DICT_TO_NUMBER = {
    'p': -21,
    'q': -22,
    'l': -23,
    'k': -24,
    'm': -25,
    'n': -26,
    'o': -27,
    'e': -28,
    'h': -29,
    'g': -10,
    'i': -11,
    'j': -12,
    'b': -13,
    'a': -14,
    'c': -15,
    'd': -16,
    'f': -17,
    '1': 1,
    '2': 2,
    '3': 3,
    '6': 6,
    '7': 7,
    '8': 8,
    '4': 4,
    '5': 5, 
    '0': None
}

"""
Translate from bonsai number ids to toys. Follow the bonsai format (positive numbers
for genotyped individuals, negative for ungenotyped)
"""
TOY_DICT_TO_STRING = {
    -21: 'p',
    -22: 'q',
    -23: 'l',
    -24: 'k',
    -25: 'm',
    -26: 'n',
    -27: 'o',
    -28: 'e',
    -29: 'h',
    -10: 'g',
    -11: 'i',
    -12: 'j',
    -13: 'b',
    -14: 'a',
    -15: 'c',
    -16: 'd',
    -17: 'f',
    1: '1',
    2: '2',
    3: '3',
    6: '6',
    7: '7',
    8: '8',
    4: '4',
    5: '5',
}

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

def translate_ids_to_nums_pedigrees(pedigree_dict):
    """
    up_dict1 = {
        1 : [None, None, -1, -2],
        2 : [None, None, -1, -2],
        3 : [None, None, -3, -4],
        4 : [None, None, -3, -4],
        -1 : [None, None, -5, -6],
        -3 : [None, None, -5, -7],
    }
    """
    new_dict = {}
    for key, value in pedigree_dict.items():
        new_key = TOY_DICT_TO_NUMBER[key]
        new_value = value
        new_value[-2] = TOY_DICT_TO_NUMBER[value[-2]]
        new_value[-1] = TOY_DICT_TO_NUMBER[value[-1]]

        if not new_value[-1] or not new_value[-2]:
            continue
        new_dict[new_key] = new_value
    
    return new_dict

def translate_ids_to_nums_ibds(ibd_lists):
    for line in ibd_lists:
        line[0] = TOY_DICT_TO_NUMBER[line[0]]
        line[1] = TOY_DICT_TO_NUMBER[line[1]]

    return ibd_lists

test_dict = translate_ids_to_nums_pedigrees(toy_dict1)

def main():
    test_list1 = translate_ids_to_nums_ibds(ibd_seg_list_truth)
    print(test_list1)

if __name__ == '__main__':
    main()