import infrared as ir
import infrared.rna as rna
import RNA
from tqdm import tqdm
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import random
import math

import argparse

parser = argparse.ArgumentParser(description="A script to create leader for multiple designs")
parser.add_argument("-s", "--steps", type=int, help="steps for MC optimization", default=100000)

args = parser.parse_args()

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

def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)

# convert extended iupac nts to RNA nts
def iupacvalues(symbol):
    return [ rna.nucleotide_to_value(x) for x in extended_iupac_nucleotides[symbol] ]

# convert IR sample to RNA sequence
def sample_to_seq(sample):
    return rna.values_to_seq(sample.values())

# calculate GC content of seq
def calculate_gc(seq):
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

# calculate target frequency of given sequence and ss
def target_frequency(sequence, ss):
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)

# check if secondary structure is unpaired in the beginning
def acceptable_seq(ss):
    if not ss.startswith('........'):
        return False
    return True

# test suboptimal secondary structures and check if they interfere with xrRNA structure or unpaired region of the leader
def test_suboptimal_ss(leader_seq, xrRNA_seq, xrRNA_ss):
    target_ss = remove_positioned_gaps(leader_seq, target)
    target_ss = target_ss + xrRNA_ss

    seq = leader_seq.replace('-','') + xrRNA_seq

    fc = RNA.fold_compound(seq)
    (ss, mfe) = fc.mfe()
    strikes = 0
    for s in fc.subopt(500):
        if s.structure[len(leader_seq):] != xrRNA_ss or not s.structure.startswith('.....'):
            return False
    return True

# print suboptimal secondary structures
def print_suboptimal_ss(leader_seq, xrRNA_seq, xrRNA_ss, sample_size=1000):
    ss_leader = remove_positioned_gaps(leader_seq, target)
    target_ss = ss_leader + xrRNA_ss

    seq = leader_seq.replace('-','') + xrRNA_seq
    
    fc = RNA.fold_compound(seq)
    (ss, mfe) = fc.mfe()
    suboptimal_ss = []
    for s in fc.subopt(sample_size):
        if s.structure not in suboptimal_ss:
            print(seq)
            print(f"{s.structure}\t{s.energy:6.2f}")
            print(target_ss)
            print(f'{"-" * len(ss_leader) + "+" * len(xrRNA_seq)}\n')
            suboptimal_ss.append(s.structure)
    return True

# test wheter there are interactions between the xrRNA and leader structure
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

# sequences from new designs
design2 = 'UGUCAGGCCUAG-GAAAAGA-CUAGCCACGGAU--------GCGUUC--A-----AUUCGUGCAGCCUGUUU-----GAGUGUUU'
design3 = 'UGUCAGGCCUCG-UCCAGU--CGAGCCACGUGCC-------GGUCCG--A----GGUACGUGCAGCCUGUUUUUC--CGGAUUUU'
design7 = 'UGUCAGGCCCUC-GUUCA---GAGGCCACGUCU-A------GAG-----CAAAA-AGACGUGCAGCCUGUU---------UUUUU'
design10 =  'UGUCAGGCCUCUGUAAC---CAGAGCCACGUGU--------UUUAACAUU-----GCACGUGCAGCCUGUUUUU-AUGUUAAGUC'

designs = [design2, design3,design7, design10]
d_nums = [2, 3, 7, 10]
structures = ['...((((((((((.......))))).((((((((....................))))))))..)))))................',
		      '((............................................................)).....................',
		      '.........................................((((((((..........................))))))))..']

xrRNA_strucs = []
xrRNA_seqs = []
for seq in designs:
	xrRNA_strucs.append([remove_positioned_gaps(seq, s) for s in structures])
	xrRNA_seqs.append(seq.replace('-',''))


## GCUAA stays unpaired
target =         ".........((((((....))))))...."
iupac_sequence = 'GCUAANNNXNNNNNXNNNXXNNNNNNNXX'
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

    for i in range(steps):  # tqdm
        cc = random.choices(ccs, weights)[0]
        new = sampler.resample(cc, cur)

        new_seq = rna.values_to_seq(new.values())
        newval = objective(new_seq)
        if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
            cur, curval = new, newval
            if curval > bestval:
                best, bestval = cur, curval
        
        print(f'Done: {i/steps*100:.2f}%\tbest TF: {bestval:.2f}', end='\r')

    return (best, bestval), sampler

def objective_function(sequence, target, freq_threshold=0.1):
    global xrRNA_strucs
    global xrRNA_seqs
    cur_sum_freq = 0
    for j in range(len(designs)):
        target_ss = remove_positioned_gaps(sequence, target)
        target_ss += xrRNA_strucs[j][0]
    
        whole_seq = sequence + xrRNA_seqs[j]
        cur_freq = target_frequency(whole_seq, target_ss)
        if cur_freq <= freq_threshold:
            return 0
        cur_sum_freq += cur_freq
    return cur_sum_freq


(best, bestval), sampler = mc_optimize(model,
                                lambda sequence: objective_function(sequence, target),
                                args.steps,
                                0.01,
                                )

print('-'*100)
# get sequence from the MC optimization
best_seq = rna.values_to_seq(best.values())
# get structure from sample
culled_structures = remove_positioned_gaps(best_seq, target)
print('best leader sequence:\t', best_seq)
culled_seq = best_seq.replace('-','')
print('culled seq + ss')
print(culled_structures)
print(culled_seq)

tfs = []
for j in range(len(designs)):
    target_ss = remove_positioned_gaps(best_seq, target)
    target_ss += xrRNA_strucs[j][0]
    whole_seq = best_seq + xrRNA_seqs[j]
    cur_tf = target_frequency(whole_seq, target_ss)
    tfs.append(cur_tf)
    print(f'Design {d_nums[j]}:{cur_tf:.2f}')
print(f'sum tf: {np.sum(tfs):.3f}\tmean tf: {np.mean(tfs):.3f}')

print('-'*100)


print('sequences and their secondary structures:')
for i in range(len(designs)):
    print(f'design {d_nums[i]}')
    target_ss = remove_positioned_gaps(best_seq, target)
    seq = culled_seq.replace('-', '') + xrRNA_seqs[i]
    print(seq)
    print(target_ss + xrRNA_strucs[i][0])
    for j in range(1,3): print(len(target_ss) * '.' + xrRNA_strucs[i][j])

