import matplotlib.pyplot as plt
# import pymol
# from pymol import cmd
import os
from itertools import combinations
import re
import pandas as pd
import numpy as np

def plot_clust_vs_all(csv_file, clustx):
    '''
    create bar plot between one cluster (clustx) and all other cluster from a given csv file --> (created with calculate_RMSD)
    '''
    df = pd.read_csv(csv_file)
    cluster_w_clustx = df[df['0'].str.contains(clustx)]
    scatter_point = {}
    for index, value in cluster_w_clustx.iterrows():
        (clust_1, clust_2) = value[1][1:-1].replace("'", "").split(", ")
        rmsd = float(value[2])
        if clust_2 != clustx:
            id = int(clust_2[-2:])
            scatter_point[id] = rmsd
        else:
            id = int(clust_1[-2:])
            scatter_point[id] = rmsd

    ids = list(scatter_point.keys())
    values = list(scatter_point.values())

    plt.bar(ids, values)
    plt.hlines(np.mean(values), 0, max(ids)+1, "red")
    plt.xlabel('cluster ID')
    plt.ylabel('RMSD')
    plt.title(f'RMSD plot of cluster {clustx[-2:]} vs all other cluster')
    plt.show()


def plot_histogram(csv_file):
    '''
    create histogramm plot of the relative frequency of the RMSD values rounded of a csv file --> (created with calculate_RMSD)
    '''
    mode = csv_file[-8:-4].upper()

    df = pd.read_csv(csv_file)
    values = df['1']
    if mode == 'RMSD':
        plt.hist(values, density=True, bins=50)
    else:
        plt.hist(values, density=True)

    plt.xlabel(mode)
    plt.ylabel('relative frequency')
    #plt.xticks([i for i in range(int(min(relative_freq.index)), int(max(relative_freq.index)) +1)])
    plt.title(f'Histogram of {mode}')
    plt.show()
    #plt.savefig(f'hist_{mode}.png')


def calculate_RMSD(qrna_folder = ".", output_csv = "pymol_rmsd.csv"):
    if not os.path.isdir('data'):
        os.system("mkdir data")
        os.system("ln -s /usr/lib64/python3.11/site-packages/pymol/data/pymol/ ./data/pymol")

    qrna_files = os.listdir(qrna_folder)
    qrna_files = [qrna_folder + "/" + file for file in qrna_files if file.endswith(".pdb")]
    mol_names = []

    for file in qrna_files:
        name = re.search("(clust\d+)_qrna", file)[1]
        mol_names.append(name)
        # print(f"cmd.load({file}, '{name}')")
        cmd.load(file, name)

    mol_combs = combinations(mol_names, 2)

    output_pd = []
    for mol1, mol2 in list(mol_combs):
        cmd.align(mol1, mol2)
        rmsd = cmd.rms_cur(mol1, mol2)
        output_pd.append([(mol1, mol2), rmsd])

    df = pd.DataFrame(output_pd)
    df.to_csv(output_csv)
    print(f"calculated RMSD were save to {output_csv}")



def calculate_SASA(qrna_folder = ".", output_csv = "pymol_sasa.csv"):
    if not os.path.isdir('data'):
        os.system("mkdir data")
        os.system("ln -s /usr/lib64/python3.11/site-packages/pymol/data/pymol/ ./data/pymol")

    qrna_files = os.listdir(qrna_folder)
    qrna_files = [qrna_folder + "/" + file for file in qrna_files if file.endswith(".pdb")]
    mol_names = []

    for file in qrna_files:
        name = re.search("(clust\d+)_qrna", file)[1]
        mol_names.append(name)
        # print(f"cmd.load({file}, '{name}')")
        cmd.load(file, name)

    output_pd = []
    for mol in mol_names:
        sasa = cmd.get_area(selection=mol)
        output_pd.append([mol, sasa])

    df = pd.DataFrame(output_pd)
    df.to_csv(output_csv)
    print(f"calculated SASA were save to {output_csv}")
    print(f"Mean SASA: {df['1'].mean()}")





min_clust = 1
max_clust = 41


def main():
    #calculate_RMSD()
    #calculate_SASA
    #histogramm("seq3_mc_800000/pymol_rmsd.csv")
    #plot_clust_vs_all("seq3_mc_800000/pymol_rmsd.csv", "clust02")

    sasa_csv = '/Users/katringutenbrunner/Desktop/MA/working/simRNA/03_05/pymol_sasa.csv'
    rmsd_csv = '/Users/katringutenbrunner/Desktop/MA/working/simRNA/03_05/pymol_rmsd.csv'
    df = pd.read_csv(sasa_csv)
    print(df["1"].mean())
    plot_histogram(rmsd_csv)
    plot_histogram(sasa_csv)



if __name__ == "__main__":
    main()