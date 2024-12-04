from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
import random
import math
import ir_utils as ir_ut
import matplotlib.pyplot as plt
import numpy as np

settings = RNA.md()
RNA.read_parameter_file('data/rna_andronescu2007.par')
#NA.read_parameter_file('data/rna_turner2004.par')


def margin_left(left_text, right_text, padding):
    # Prints two text (left and right) with a certain padding between them
    print(f"{left_text: <{padding}}{right_text}")


def sigmoid(x):
    # Calculates the sigmoid function of input
    return 1 / (1 + math.exp(-x))

def remove_positioned_gaps(sequence, structure):
    # Removes the gaps of given structures based on the given sequence
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def add_gaps(sequence, structure):
    # Adds the gaps back into given structures based on the given sequence
    structure = list(structure)
    new_ss = [structure.pop(0) if nt != '-' else '-' for  nt in sequence]
    return ''.join(new_ss)


def target_frequency(sequence, model_input):
    # calculates the target frequency of a given sequence based on the structure of the given ModelInput
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)


def ensemble_defect(sequence, model_input):
    # calculates the ensemble defect of a given sequence based on the structure of the given ModelInput
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return - fc.ensemble_defect(ss)


def print_pk2_energy(sequence, model_input):
    # prints the energy of pseudoknot 2 given a sequence based on the structure of the given ModelInput
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    ss_pk2 = remove_positioned_gaps(sequence, model_input.structures[2])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    pk2_e = fc.eval_structure(ss_pk2)    
    print(pk2_e)
    return 


def mc_optimize(model, model_input, objective, steps, temp, start=None, verbose=False):
    """
    Args:
        model (ir.infrared.Model):
            The Infrared model object defining the constraints and structure of the optimization.
        model_input (ir_utils.ModelInput):
            Input data for the model, including structures and IUPAC sequence.
        objective (callable): 
            Objective function that evaluates the quality of a sequence
        steps (int):
            The number of Monte Carlo optimization steps to perform.
        temp (float): 
            Temperature for the MC Optimization
        start (ir.libinfrared.Assignment): 
            The initial state for the optimization. If None, a random sample is used as the 
            starting point. (default: None)
        verbose (bool):
            Dis- or enables output of optimization-process (default: False)

    Returns:
        tuple: A tuple containing:
            - best_result (tuple): A tuple `(best, bestval)` where:
                - best (ir.libinfrared.Assignment): The best state found during optimization.
                - bestval (float): The objective value of the best state.
            - sampler (ir.infrared.MultiDimensionalBoltzmannSampler): The sampler object used for the optimization.
    """

    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start
    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    # 9 is first gap in component, simple but just for now i will keep it like that
    #weights = [1 if not 9 in cc else 1/10 for cc in ccs]
    weights = [1/len(cc) for cc in ccs]
    if verbose:
        for i in tqdm(range(steps)):  # tqdm
            cc = random.choices(ccs, weights)[0]
            new = sampler.resample(cc, cur)

            new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
            newval = objective(new_seq)
            
            if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
                cur, curval = new, newval
                if curval > bestval:
                    
                    best, bestval = cur, curval
                    best_seq = rna.values_to_seq(best.values()[:len(model_input.structures[0])])
    else:
        for i in range(steps):
            cc = random.choices(ccs, weights)[0]
            new = sampler.resample(cc, cur)

            new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
            newval = objective(new_seq)
            
            if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
                cur, curval = new, newval
                if curval > bestval:
                    
                    best, bestval = cur, curval
                    best_seq = rna.values_to_seq(best.values()[:len(model_input.structures[0])])

    return (best, bestval), sampler

