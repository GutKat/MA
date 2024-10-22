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


# def neighbor_optimize_top10(sampler, model, model_input, objective, steps, temp, shift, variance, dist_func, start=None):
#     n = len(model_input.structures[0])
#     curSampler = sampler
#     top_n = 10
#     steps = 50

#     if start == None:
#         curSample = curSampler.targeted_sample()
#         cur_seq = rna.values_to_seq(curSample.values()[:n])
#         curval = objective(cur_seq)

#         best, bestval = curSample, curval
#         oldVals = curSample.values()[:n]
#         curSampler = ir.getNeighborSampler(curSampler, oldVals, dist_func, shift, variance)

#         top_neighbors = [best]
#         for _ in range(top_n-1):
#             new_neighbor = curSampler.targeted_sample()
#             top_neighbors.append(new_neighbor)

#     else:
#         curSample = start[0]
#         cur_seq = rna.values_to_seq(curSample.values()[:n])

#         best, bestval = curSample, curval
#         oldVals = curSample.values()[:n]

#         top_neighbors = start

#     old_top = top_neighbors
           

#     for _ in range(steps):
#         neighbors = [*top_neighbors]

#         for neighbor in top_neighbors:
#             neighborVal = neighbor.values()[:n]
#             ir.updateNeighborSampler(curSampler, neighborVal, oldVals, dist_func, shift, variance)
#             oldVals = neighborVal

#             for _ in range(top_n):
#                 new_neighbor = curSampler.targeted_sample()
#                 neighbors.append(new_neighbor)

#         heap = []
#         for neighbor in neighbors:
#             neighbor_seq = rna.values_to_seq(neighbor.values()[:n])
#             neighbor_val = objective(neighbor_seq)
            
#             if len(heap) < top_n:
#                 heapq.heappush(heap, (neighbor_val, neighbor))
#             else:
#                 heapq.heappushpop(heap, (neighbor_val, neighbor))

        
#         top_elements = heapq.nlargest(top_n, heap)
#         top_neighbors = [item for _, item in top_elements]
#         best, best_seq, bestval = top_elements[0][1], rna.values_to_seq(top_elements[0][1].values()[:n]), top_elements[0][0]

#         print(redundant_elements(old_top, top_neighbors))
#         old_top = top_neighbors

#         print(best_seq, len(best_seq.replace('-','')), round(bestval,4))
    
#     sampler = curSampler
#     return (best, bestval), top_neighbors, sampler


def sort_tuple(tup): 
 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using second element of 
    # sublist lambda has been used 
    return sorted(tup, key = lambda x: x[1], reverse=True)


def sort_tuple(tup): 
 
    # reverse = None (Sorts in Ascending order) 
    # key is set to sort using second element of 
    # sublist lambda has been used 
    return sorted(tup, key = lambda x: x[1], reverse=True)



def neighbor_optimize_top10(sampler, model, model_input, objective, steps, temp, shift, variance, dist_func, start=None):
    n = len(model_input.structures[0])
    curSampler = sampler
    top_n = 10
    steps = 25
    if start == None:
        curSample = curSampler.targeted_sample()
        cur_seq = rna.values_to_seq(curSample.values()[:n])
        oldVals = curSample.values()[:n]
        curSampler = ir.getNeighborSampler(curSampler, oldVals, dist_func, shift, variance)

        top_neighbors = [curSample]
        for _ in range(top_n-1):
            new_neighbor = curSampler.targeted_sample()
            top_neighbors.append(new_neighbor)
        
        neighbors = [(neigh, objective(rna.values_to_seq(neigh.values()[:n]))) for neigh in top_neighbors]
        top_neighbors = neighbors


    else:
        curSample = start[0][0]
        oldVals = curSample.values()[:n]
        top_neighbors = neighbors = start

    # old_top = top_neighbors
           

    for _ in range(steps):
        neighbors = top_neighbors[:int(len(top_neighbors))]


        for neighbor in top_neighbors:
            neighborVal = neighbor[0].values()[:n]
            ir.updateNeighborSampler(curSampler, neighborVal, oldVals, dist_func, shift, variance)
            oldVals = neighborVal

            for _ in range(top_n):
                new_neighbor = curSampler.targeted_sample()

                if len(neighbors) < 10:
                    neigh_seq = rna.values_to_seq(new_neighbor.values()[:n])
                    neighbors.append((new_neighbor, objective(neigh_seq)))
                    neighbors = sort_tuple(neighbors)

                else:
                    neigh_seq = rna.values_to_seq(new_neighbor.values()[:n])
                    neigh_obj = objective(neigh_seq) 
                    if neigh_obj > neighbors[-1][1]:
                        neighbors[-1] = (new_neighbor, neigh_obj)
                        neighbors = sort_tuple(neighbors)
        
        
        top_neighbors = neighbors
        
        
        best, bestval = top_neighbors[0]
        best_seq = rna.values_to_seq(best.values()[:n])
        print(best_seq, len(best_seq.replace('-','')), round(bestval,4))
    
    sampler = curSampler
    return (best, bestval), top_neighbors, sampler



def redundant_elements(list1, list2):
    # Convert both lists to sets to remove duplicates and then find intersection
    set1 = set(list1)
    set2 = set(list2)
    
    # Find the common elements
    common_elements = set1.intersection(set2)
    
    # Return the number of common elements
    return len(common_elements)


