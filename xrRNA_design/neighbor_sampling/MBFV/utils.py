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


def ensemble_defect(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return - fc.ensemble_defect(ss)


def neighbor_optimize(sampler, objective, steps, temp, shift, variance, dist_func, start=None):
    """Optimize by Monte-Carlo optimization with neighbor sampling

    Maximizes an objective function over assignments of a Infrared model
    using a Monte-Carlo stochastic optimization strategy with
    neighbor resampling and Metropolis criterion.
    Targeted sampling is used.

    In each iteration, the algorithm samples for assignments which are in the
    neighborhood of the currently selected assignment.
    The distance function used to define the neighbohood as well as the extent of the
    neighborhood must be provided. For RNA design, the hamming distance would be arecommendation

    If no start assignment is given, a new one will be sampled from the given sampler.

    Args:
        sampler: Infrared sampler describing assignments
        objective: objective function on assignments
        steps: iterations
        temp: temperature for the Metropolis criterion
        shift: The shift allowed for the neighborhood sampling. Given as a fraction of the number of variables.
        variance: The variance allowed for the neighborhood sampling. Given as a fraction of the number of variables.
        distanceFunction: A distance funtion which will be applied to each variable. Must take the variable of the assignemnt it compares to as an additional argument.
        start: optional start assignment

    Returns:
        Pair of best assignment and its objective value
    """

    currSampler = sampler
    currSample = currSampler.targeted_sample() if start is None else start

    cureval = objective(currSample)
    best, besteval = currSample, cureval
    newSample, neweval = currSample, cureval
    oldVals = currSample.values()

    currSampler = getNeighborSampler(currSampler, currSample.values(), dist_func, shift, variance)
    currSampler.sample()
    for _ in range(steps):
        if (currSample == newSample):
            updateNeighborSampler(currSampler, currSample.values(), oldVals, dist_func, shift, variance)
            oldVals = currSample.values()
        newSample = currSampler.targeted_sample()
        neweval = objective(newSample)
        if (neweval >= cureval or random.random() <= math.exp((neweval-cureval)/temp)):
            currSample, cureval = newSample, neweval
            if cureval > besteval:
                best, besteval = currSample, cureval

    return (best, besteval)