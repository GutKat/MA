import MDAnalysis
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import math
import os
import sys
import statistics
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from pylab import *

def get_pairwise_distances(dcd_file, pdb_file, movie_skip):
    u = MDAnalysis.Universe(pdb_file, dcd_file)
    selection = u.select_atoms("record_type ATOM")
    pairwise_distances = []

    for t, ts in enumerate(u.trajectory):
        if t % movie_skip != 0:
            continuev.add_aux_BP(14, 20, edge5="s", color="#FF00FF")
        positions = selection.center(None, compound="residues")

        for i, pos_i in enumerate(positions):
            pairwise_distances[t].append([])
            for j, pos_j in enumerate(positions):
                if i > j:
                    pairwise_distances[t][i].append(pairwise_distances[t][j][i])
                else:
                    pairwise_distances[t][i].append(np.linalg.norm(pos_j-pos_i))

    return pairwise_distances

def get_bp_distance(dcd_file, pdb_file, movie_skip, bp):
    u = MDAnalysis.Universe(pdb_file, dcd_file)
    selection = u.select_atoms("record_type ATOM")
    distances = []

    for t, ts in enumerate(u.trajectory):
        if t % movie_skip != 0:
            continue
        t = int(t/movie_skip)
        positions = selection.center(None, compound="residues")
        pos_i = positions[bp[0]]
        pos_j = positions[bp[1]]
        distances.append(np.linalg.norm(pos_j-pos_i))

    return distances

def remove_comments(lines):
    cleaned = [line for line in lines if line[0].isnumeric()]
    return cleaned

def get_steps(foldername, gradient):
    outLines = []
    inF = open(f"{foldername}/out{gradient}_1_1.txt", 'r')

    outLines.append("Step ")
    lines = remove_comments(inF.readlines())

    start = int(lines[0][0:lines[0].find(',')])
    step = int(lines[1][0:lines[1].find(',')]) - int(lines[0][0:lines[0].find(',')])
    inF.close()

    return start, step


def calculate_endforce(foldername, distances, ):
    end_force = ((len(distances)*movie_skip*step_interval*2)/1000000000) * gradient
    return end_force

base_folder = '/scr/aldea/kgutenbrunner/working/MD/MBFV/designs/design_02/'
foldername = os.path.join(base_folder, 'run1000_1/')
gradient = 1000

movie_skip = 50

basedir = os.getcwd()

dcd_files = [foldername + i for i in os.listdir(foldername) if i.endswith('.dcd')]
pdb_file = base_folder + 'init.pdb'



step_start, step_interval = get_steps(foldername, gradient)

distances_file = os.path.join(base_folder, 'distances.txt')

distances_5_end = {}
with open(distances_file, 'r') as f:
    distances = None
    simulation = None
    for line in f.readlines():
        if not line.startswith('>'):
            distances.append(float(line))
        else:
            distances_5_end[simulation] = distances
            simulation=int(line[1:-1]) +1
            distances=[]

del distances_5_end[None]
distances_5_end[simulation] = distances
force=line

for key,distances in distances_pk1.items():
    end_force = calculate_endforce(foldername, distances)
    plt.plot(np.linspace(start=0, stop=end_force, num=len(distances[1:])), distances[1:], label=key)

plt.legend(title='Simulation', loc=4)
plt.title(f'5\' End Distance during MD simulations', fontsize=14)
plt.xlabel('Force [pN]')
plt.ylabel(r'Distance [$\AA$]')

save_under = os.path.join(base_folder, 'distance_5end.png')
print(save_under)
plt.savefig(save_under)
