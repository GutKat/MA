from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
from collections import namedtuple
import random
import math
import ir_utils as ir_ut
import matplotlib.pyplot as plt
import numpy as np


def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")


target_len = 87
target_gc =  0.58
target_energy = -30
energy_range = [-35, -20]
relations_median= {
        'st2_hl2': 3,
        'st3_hl3': 2,
        'sI_sII': 1.823529,
        'sII_sIII': 0.705882,
        'sI_sIII': 1.35
    }

relations_range = {
        'st2_hl2': [2, 4],
        'st3_hl3': [3/4, 2],
        'sI_sII': [1.5, 2.6],
        'sII_sIII': [0.5, 0.9],
        'sI_sIII': [1.20, 1.4]
    }

penalty = 0.3
length_fraction = 0.0
relation_fraction = 0.0 
plotting_steps = 100


def number_in_range(x, range_):
    if range_[0] <= round(x, 3) <= range_[1]:
        return True
    return False


def error(x, y):
    return abs(x-y)


def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def add_gaps(sequence, structure):
    structure = list(structure)
    new_ss = [structure.pop(0) if nt != '-' else '-' for  nt in sequence]
    return ''.join(new_ss)


def target_frequency(sequence, model_input):
    target = model_input.structures[0]
    ss = remove_positioned_gaps(sequence, target)
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)


def mc_optimize(model, model_input, objective, steps, temp, start=None):
    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start
    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1/len(cc) for cc in ccs]




    for i in ((range(steps))): #tqdm
        cc = random.choices(ccs,weights)[0]
        new = sampler.resample(cc, cur)
        new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
        newval = objective(new_seq)

    
        if (newval >= curval or random.random() <= math.exp((newval-curval)/temp)):
            cur, curval = new, newval
            if curval > bestval:
                best, bestval = cur, curval


    return (best, bestval), sampler

