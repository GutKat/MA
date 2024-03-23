from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
from collections import namedtuple
import random
import math



def target_frequency(sequence, target):
    ss = remove_positioned_gaps(sequence, target)
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)


def remove_positioned_gaps(sequence, structure):
    print(sequence)
    print(structure)
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)

    


def mc_optimize(model, model_input, objective, steps, temp, start=None):
    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start

    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1/len(cc) for cc in ccs]
    for i in tqdm(range(steps)):
        cc = random.choices(ccs,weights)[0]
        new = sampler.resample(cc, cur)
        new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
        newval = objective(new_seq)
        if (newval >= curval or random.random() <= math.exp((newval-curval)/temp)):
            cur, curval = new, newval
            if curval > bestval:
                #print(curval)
                best, bestval = cur, curval
    return (best, bestval)



