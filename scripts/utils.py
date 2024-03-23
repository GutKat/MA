import numpy as np


def extract_db_structure(ss_file):
    with open(ss_file, 'r') as f:
        true_structure = f.read()
    
    true_structure = true_structure.split('\n')
    db_dict = {0:"()", 1:'[]', 2:'{}', 3:'<>'}
    db_ss = list('.' * len(true_structure[0]))

    for i, ss in enumerate(true_structure):
        for j in range(len(ss)):
            if ss[j] == '(':
                db_ss[j] = db_dict[i][0]
            if ss[j] == ')':
                db_ss[j] = db_dict[i][1]

    return ''.join(db_ss)




def db_to_matrix(db):
    matrix = np.zeros((len(db), len(db)), dtype=int)
    db_dict = {')': '(',
               ']': '[',
               '}': '{',
               '>': '<'}
    db_stacks = {'(':[],
               '[':[],
               '{':[],
               '<':[],
               }

    for i, symbol in enumerate(db):
        if symbol in db_stacks.keys():
            db_stacks[symbol].append(i)
        elif symbol in ')}]>':
            opening_index = db_stacks[db_dict[symbol]].pop()
            closing_index = i
            matrix[opening_index][closing_index] = 1
            matrix[closing_index][opening_index] = 1
    return matrix

