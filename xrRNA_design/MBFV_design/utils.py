from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
import random
import math
import ir_utils as ir_ut
import matplotlib.pyplot as plt
import numpy as np


def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")


def sigmoid(x):
    return 1 / (1 + math.exp(-x))

def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def add_gaps(sequence, structure):
    structure = list(structure)
    new_ss = [structure.pop(0) if nt != '-' else '-' for  nt in sequence]
    return ''.join(new_ss)


def target_frequency(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)




def objective_function(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    ss_pk2 = remove_positioned_gaps(sequence, model_input.structures[2])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    ensemble_defect = fc.ensemble_defect(ss)
    freq = fc.pr_structure(ss)

    pk2_e = fc.eval_structure(ss_pk2)    
    pk2_e_sig = sigmoid(pk2_e)
    result = freq - (ensemble_defect + pk2_e_sig)
    return result


def mc_optimize(model, model_input, objective, steps, temp, start=None):
    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start
    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1 / len(cc) for cc in ccs]
    for i in tqdm((range(steps))):  # tqdm
        cc = random.choices(ccs, weights)[0]
        new = sampler.resample(cc, cur)

        new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
        newval = objective(new_seq)
        if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
            cur, curval = new, newval
            if curval > bestval:
                best, bestval = cur, curval

    return (best, bestval), sampler


######################################
# for testing the objective function #
######################################

def objective_function0(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    ss_pk2 = remove_positioned_gaps(sequence, model_input.structures[2])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    freq = fc.pr_structure(ss)

    pk2_e = fc.eval_structure(ss_pk2)    
    pk2_e_sig = sigmoid(pk2_e)
    result = freq - pk2_e_sig
    return result

def objective_function1(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    ss_pk2 = remove_positioned_gaps(sequence, model_input.structures[2])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    ensemble_defect = fc.ensemble_defect(ss)

    pk2_e = fc.eval_structure(ss_pk2)    
    pk2_e_sig = sigmoid(pk2_e)
    result = - (ensemble_defect + pk2_e_sig)
    return result


def objective_function2(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    ss_pk2 = remove_positioned_gaps(sequence, model_input.structures[2])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    ensemble_defect = fc.ensemble_defect(ss)
    freq = fc.pr_structure(ss)

    pk2_e = fc.eval_structure(ss_pk2)    
    pk2_e_sig = sigmoid(pk2_e)
    result = freq - (ensemble_defect + pk2_e_sig)
    return result
