import infrared as ir
import infrared.rna as rna
import RNA
from tqdm import tqdm
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import random
import math

def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)

extended_iupac_nucleotides = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'U': 'U',
    'R': 'AG',
    'N': 'ACGU',
    'X': 'ACGU-',
    '-': '-',
    '.': '-'
}

def sample_to_seq(sample):
    return rna.values_to_seq(sample.values())

def calculate_gc(seq):
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

def iupacvalues(symbol):
    return [ rna.nucleotide_to_value(x) for x in extended_iupac_nucleotides[symbol] ]

def test_interaction(leader_seq, xrRNA_seq, xrRNA_ss,threshold = 0.01):
    leader_seq = leader_seq.replace('-','')

    seq = leader_seq + xrRNA_seq
    fold_compound = RNA.fold_compound(seq)
    fold_compound.pf()
    base_pair_probs = np.array(fold_compound.bpp())
    for i in range(len(leader_seq)):
        for j in range(len(leader_seq), len(seq)):
            if base_pair_probs[i][j] > threshold:  # Adjust threshold as needed
                return False
    return True


def target_frequency(sequence, ss):
    ss = remove_positioned_gaps(sequence, ss)
    whole_seq = sequence + xrRNA_seq
    whole_ss = ss + xrRNA_ss[0]
    fc = RNA.fold_compound(whole_seq.replace('-',''))
    fc.pf()
    return fc.pr_structure(whole_ss)

def print_suboptimal_ss(leader_seq, xrRNA_seq, xrRNA_ss, sample_size=1000):
    target_ss = remove_positioned_gaps(leader_seq, target)
    
    seq = leader_seq.replace('-','')
    seq = seq + xrRNA_seq
    
    target_ss = target_ss + xrRNA_ss
    
    fc = RNA.fold_compound(seq)
    (ss, mfe) = fc.mfe()
    suboptimal_ss = []
    for s in fc.subopt(sample_size):
        if s.structure not in suboptimal_ss:
            print(seq)
            print(f"{s.structure}\t{s.energy:6.2f}")
            print(target_ss)
            print(f'{"-" * len(leader_seq) + "+" * len(xrRNA_seq)}\n')
            suboptimal_ss.append(s.structure)
    return True

seq = 'AGUCAGGCCGGG-UCC-----CCCGCCACGUGGAG------CCCU----UA---UCCGCGUGCUGCCUGU---------AGGGAA'

structures = ['...((((((((((.......))))).((((((((....................))))))))..)))))................',
		      '((............................................................)).....................',
		      '.........................................((((((((..........................))))))))..']
xrRNA_ss = [remove_positioned_gaps(seq, s) for s in structures]
xrRNA_seq = seq.replace('-','')

## GCUAA stays unpaired
target =         ".........(((((.....)))))...."
iupac_sequence = 'GCUAANNNXNNNNNNNNXXNNNNNNNXX'
n = len(target)
model = ir.Model()
model.add_variables(n, 5)

for i, x in enumerate(iupac_sequence):
    model.add_constraints(ir.ValueIn(i, iupacvalues(x)))

model.add_constraints(rna.BPComp(i,j) for (i,j) in rna.parse(target))
model.add_functions([rna.GCCont(i) for i in range(n)], 'gc')



def mc_optimize(model, objective, steps, temp):
    sampler = ir.Sampler(model)
    cur = sampler.sample()
    cur_seq = rna.values_to_seq(cur.values())
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1 / len(cc) for cc in ccs]

    for i in tqdm((range(steps))):  # tqdm
        cc = random.choices(ccs, weights)[0]
        new = sampler.resample(cc, cur)

        new_seq = rna.values_to_seq(new.values())
        newval = objective(new_seq)
        if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
            cur, curval = new, newval
            if curval > bestval:
                best, bestval = cur, curval
    return (best, bestval), sampler

steps = 100000
(best, bestval), sampler = mc_optimize(model,
                                lambda sequence: target_frequency(sequence, target),
                                steps,
                                0.01,
                                )

print(bestval)

# get sequence from the MC optimization
best_seq = rna.values_to_seq(best.values())
# get structure from sample
culled_structures = remove_positioned_gaps(best_seq, target)
culled_seq = best_seq.replace('-','')
print(culled_structures)
print(culled_seq)

print('\n')
print('-'*100)
print('\n')

print(best_seq)
print_suboptimal_ss(best_seq, xrRNA_seq, xrRNA_ss[0], 400)
print('\n')
print('-'*100)
print('\n')
