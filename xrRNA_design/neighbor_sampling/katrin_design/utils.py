from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
import random
import math
import ir_utils as ir_ut
import matplotlib.pyplot as plt
import numpy as np
import heapq


def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")


def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def add_gaps(sequence, structure):
    structure = list(structure)
    new_ss = [structure.pop(0) if nt != '-' else '-' for nt in sequence]
    return ''.join(new_ss)


def target_frequency(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-', ''))
    fc.pf()
    return fc.pr_structure(ss)


def ensemble_defect(sequence, model_input):
    ss = remove_positioned_gaps(sequence, model_input.structures[0])
    fc = RNA.fold_compound(sequence.replace('-', ''))
    fc.pf()
    return - fc.ensemble_defect(ss)


def neighbor_optimize_top10(sampler, model, model_input, objective, steps, temp, shift, variance, dist_func, start=None):
    n = len(model_input.structures[0])
    curSampler = sampler
    top_n = 10
    steps = 50

    if start == None:
        curSample = curSampler.targeted_sample()
        cur_seq = rna.values_to_seq(curSample.values()[:n])
        curval = objective(cur_seq)

        best, bestval = curSample, curval
        oldVals = curSample.values()[:n]
        curSampler = ir.getNeighborSampler(curSampler, oldVals, dist_func, shift, variance)

        top_neighbors = []
        for _ in range(top_n-1):
            new_neighbor = curSampler.targeted_sample()
            top_neighbors.append(new_neighbor)
        neighbors = [best, *top_neighbors] 

    else:
        curSample = start[0]
        oldVals = curSample.values()[:n]
        top_neighbors = start
        neighbors = [*top_neighbors] 

    old_top = top_neighbors
           

    for _ in range(steps):
        for neighbor in top_neighbors:
            # print(rna.values_to_seq(neighbor.values()[:len(model_input.structures[0])]))
            # neighborSampler = ir.Sampler(curSampler.model)
            # neighborSampler = curSampler.copySampler()

            neighborVal = neighbor.values()[:n]
            ir.updateNeighborSampler(curSampler, neighborVal, oldVals, dist_func, shift, variance)
            
            for _ in range(top_n):
                # new_neighbor = curSampler.targeted_sample()
                new_neighbor = curSampler.sample()
                # print(rna.values_to_seq(new_neighbor.values()[:len(model_input.structures[0])]))
                neighbors.append(new_neighbor)

            oldVals = neighborVal


        heap = []
        for neighbor in neighbors:
            neighbor_seq = rna.values_to_seq(neighbor.values()[:n])
            neighbor_val = objective(neighbor_seq)
            
            if len(heap) < top_n:
                heapq.heappush(heap, (neighbor_val, neighbor))
            else:
                heapq.heappushpop(heap, (neighbor_val, neighbor))

        
        top_elements = heapq.nlargest(top_n, heap)
        top_neighbors = [item for _, item in top_elements]
        best, best_seq, bestval = top_elements[0][1], rna.values_to_seq(top_elements[0][1].values()[:n]), top_elements[0][0]

        print(redundant_elements(old_top, top_neighbors))
        old_top = top_neighbors

        bestVals = best.values()[:n]
        ir.updateNeighborSampler(curSampler, bestVals, oldVals, dist_func, shift, variance)
        oldVals = bestVals

        print(best_seq, len(best_seq.replace('-','')), round(bestval,4))
    
    sampler = curSampler
    return (best, bestval), top_neighbors, sampler


# def neighbor_optimize_top10(sampler, model, model_input, objective, steps, temp, shift, variance, dist_func, start=None):
#     n = len(model_input.structures[0])

#     curSampler = sampler
#     if start is None:
#         top_neighbors = [curSampler.targeted_sample() for _ in range(10)] 
#     else:
#         top_neighbors = start

#     bestval = -100
#     best = None
#     for neighbor in top_neighbors:
#         neighbor_seq = rna.values_to_seq(neighbor.values()[:n])
#         neighbor_val = objective(neighbor_seq)
        
#         if neighbor_val > bestval:
#             bestval = neighbor_val
#             best = neighbor

#     # newSample, neweval = curSample, curval
#     oldVals = best.values()[:n]
#     curSampler = ir.getNeighborSampler(curSampler, oldVals, dist_func, shift, variance)

#     old_top = top_neighbors
#     neighbors = []

#     for _ in range(2):
#         for neighbor in top_neighbors:
#             # print(rna.values_to_seq(neighbor.values()[:len(model_input.structures[0])]))
#             # neighborSampler = ir.Sampler(curSampler.model)
#             # neighborSampler = curSampler.copySampler()

#             neighborVal = neighbor.values()[:n]
#             ir.updateNeighborSampler(curSampler, neighborVal, oldVals, dist_func, shift, variance)
            
#             for _ in range(10):
#                 # new_neighbor = curSampler.targeted_sample()
#                 new_neighbor = curSampler.sample()
#                 # print(rna.values_to_seq(new_neighbor.values()[:len(model_input.structures[0])]))
#                 neighbors.append(new_neighbor)

#             oldVals = neighborVal

#         heap = []
#         for neighbor in neighbors:
#             neighbor_seq = rna.values_to_seq(neighbor.values()[:n])
#             neighbor_val = objective(neighbor_seq)
            
#             if len(heap) < 10:
#                 heapq.heappush(heap, (neighbor_val, neighbor))
#             else:
#                 heapq.heappushpop(heap, (neighbor_val, neighbor))

        
#         top_elements = heapq.nlargest(10, heap)
#         top_neighbors = [item for _, item in top_elements]
#         best, best_seq, bestval = top_elements[0][1], rna.values_to_seq(top_elements[0][1].values()[:n]), top_elements[0][0]


#         print(redundant_elements(old_top, top_neighbors))
#         old_top = top_neighbors

#         # bestVals = best.values()[:n]
#         # ir.updateNeighborSampler(curSampler, bestVals, oldVals, dist_func, shift, variance)
#         # oldVals = bestVals

#         print(best_seq, len(best_seq.replace('-','')), round(bestval,4))

#     return (best, bestval),top_neighbors,  sampler
        

def redundant_elements(list1, list2):
    # Convert both lists to sets to remove duplicates and then find intersection
    set1 = set(list1)
    set2 = set(list2)
    
    # Find the common elements
    common_elements = set1.intersection(set2)
    
    # Return the number of common elements
    return len(common_elements)



def neighbor_optimize(sampler, model, model_input, objective, steps, temp, shift, variance, dist_func, start=None):
    curSampler = sampler
    curSample = curSampler.targeted_sample() if start is None else start

    cur_seq = rna.values_to_seq(curSample.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)

    best, bestval = curSample, curval
    newSample, neweval = curSample, curval
    oldVals = curSample.values()

    curSampler = ir.getNeighborSampler(curSampler, curSample.values(), dist_func, shift, variance)
    curSampler.sample()
    # for i in (range(steps)):
    for i in tqdm(range(steps)):
        if (curSample == newSample):
            ir.updateNeighborSampler(curSampler, curSample.values(), oldVals, dist_func, shift, variance)
            oldVals = curSample.values()
        newSample = curSampler.targeted_sample()
        new_seq = rna.values_to_seq(newSample.values()[:len(model_input.structures[0])])
        newval = objective(new_seq)

        if (newval >= curval or random.random() <= math.exp((newval - curval) / temp)):
            # cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
            # print(new_seq == cur_seq)
            curSample, curval = newSample, neweval
            if curval > bestval:

                best, bestval = curSample, curval
    #            print(bestval)
                best_seq = rna.values_to_seq(best.values()[:len(model_input.structures[0])])

    return (best, bestval), sampler
