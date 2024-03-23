import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utils as ut


def summary_analysis(simRNA_folder, file_name = False, plot_cluster_distribution = False):
    summary_file = simRNA_folder + '/' + 'summary.txt'

    # search for conensus file
    if file_name:
        structure_file = simRNA_folder + '/' + 'simRNA_' + file_name  + '_xrRNA' + '.ss_con'
    else:
        structure_file = simRNA_folder + '/' + 'simRNA_' + os.path.basename(simRNA_folder)  + '_xrRNA' + '.ss_con'

    # read in the summary file with all the secondary structures
    with open(summary_file, 'r') as f:
        data = f.read()
    data = data.split('\n')
    names = [re.search('clust\d*', data[i])[0] for i in range(0,len(data), 3) if data[i]]
    ss_dbs = [(data[i]) for i in range(2,len(data), 3) if data[i]]

    # convert dotbrackets of ss to matrix
    ss_matrices = [ut.db_to_matrix(ss) for ss in ss_dbs]

    # get and convert dotbracket of 'true structure' to matrix
    db_ss = ut.extract_db_structure(structure_file)
    true_structure = ut.db_to_matrix(''.join(db_ss))

    # get number of correct ss
    true_pred = [1 if np.array_equal(ss, true_structure) else 0 for ss in ss_matrices]

    # create dataframe of name, ss in dotbrackets and if ss is correct
    df = pd.DataFrame(list(zip(names, ss_dbs, true_pred)), columns=['name', 'db', 'correct'])
    print(f'Total structure:\t{len(ss_dbs)}\nFolding correct:\t{sum(true_pred)}\nin %:\t{round(sum(true_pred)/len(ss_dbs), 2)}')

    # get number of simulation in different clusters
    if plot_cluster_distribution:
        cluster_folder = [file for file in os.listdir(simRNA_folder) if file.startswith('all_thrs5.00A_clust') and os.path.isdir(file)]
        n_files = [len([file for file in os.listdir(folder) if file.endswith('_AA.pdb')]) for folder in cluster_folder]
        names = [i for i in range(1, len(cluster_folder) + 1)]
        rel_n_files = [i / sum(n_files) for i in n_files]
        plt.bar(names, rel_n_files)
        plt.xticks(names)
        plt.xlabel('Cluster')
        plt.ylabel('Number of simulations in cluster')
        plt.savefig('cluster.png')
    return df, [true_structure, ss_matrices]


def summary_analysis_virus(simRNA_folder, save = False):
    summary_file = simRNA_folder + '/' + 'summary.txt'

    with open(summary_file, 'r') as f:
        data = f.read()
    data = data.split('\n')
    ss_dbs = [(data[i]) for i in range(2, len(data), 3) if data[i]]
    ss_matrices = [ut.db_to_matrix(ss) for ss in ss_dbs]
    if save:
        bp_plot(ss_matrices, simRNA_folder)
    else:
        bp_plot(ss_matrices)
    return


def bp_plot(structs, folder = False):
    bp_matrix = np.sum(structs, axis = 0) / len(structs)

    #result_matrix = result_matrix / len(structs)
    plt.imshow(bp_matrix, cmap='viridis', interpolation='nearest')
    plt.xticks(np.arange(0, bp_matrix.shape[1], 5), np.arange(0, bp_matrix.shape[1], 5))
    plt.yticks(np.arange(0, bp_matrix.shape[0], 5), np.arange(0, bp_matrix.shape[0], 5))
    plt.title('Base pairs')
    plt.grid(color='gray', linestyle='--', linewidth=0.4)
    plt.colorbar()  # Add color bar
    if folder:
        plt.savefig(folder + '/' + 'heatmap_basepairs.png')
        plt.clf()
    else:
        plt.show()
        plt.clf()

def bp_plot_vs_ts(structs, true_ss, folder = False):
    fig, axs = plt.subplots(1, 2, figsize=(10,5))
    bp_matrix = np.sum(structs, axis = 0) / len(structs)

    im = axs[0].matshow(true_ss)
    axs[0].set_title('Base pairs of true structure')
    axs[0].set_xticks(np.arange(0, bp_matrix.shape[1], 10), np.arange(0, bp_matrix.shape[1], 10))
    axs[0].set_yticks(np.arange(0, bp_matrix.shape[0], 10), np.arange(0, bp_matrix.shape[0], 10))
    axs[0].grid(color='gray', linestyle='--', linewidth=0.4)

    im = axs[1].matshow(bp_matrix, cmap='viridis', interpolation='nearest')
    axs[1].set_xticks(np.arange(0, bp_matrix.shape[1], 10), np.arange(0, bp_matrix.shape[1], 10))
    axs[1].set_yticks(np.arange(0, bp_matrix.shape[0], 10), np.arange(0, bp_matrix.shape[0], 10))
    axs[1].set_title('Base pairs of clusters')
    axs[1].grid(color='gray', linestyle='--', linewidth=0.4)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    if folder:
        plt.savefig(folder + '/' + 'heatmap_basepairs.png')
        plt.clf()
    else:
        plt.show()
        plt.clf()



def bp_diff_plot(true_ss,structs, folder = False):
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    all_diffs = np.zeros_like(true_ss)
    incorrect_bps = np.zeros_like(true_ss)
    missing_bps = np.zeros_like(true_ss)

    for ss in structs:
        diff = (true_ss != ss).astype(int)
        all_diffs += diff

        diff = ((ss == 1) & (true_ss != 1)).astype(int)
        incorrect_bps += diff

        diff = ((ss == 0) & (true_ss == 1)).astype(int)
        missing_bps += diff

    all_diffs = all_diffs / len(structs)
    incorrect_bps = incorrect_bps / len(structs)
    missing_bps = missing_bps / len(structs) 
    rel_distance = np.sum(all_diffs)


    im = axs[0].matshow(all_diffs, cmap='viridis', interpolation='nearest', vmin=0, vmax=1) #, vmin=0, vmax=1
    axs[0].set_xticks(np.arange(0, all_diffs.shape[1], 10), np.arange(0, all_diffs.shape[1], 10))
    axs[0].set_yticks(np.arange(0, all_diffs.shape[0], 10), np.arange(0, all_diffs.shape[0], 10))
    axs[0].set_title('Total differences of clusters and true structure')
    axs[0].grid(color='gray', linestyle='--', linewidth=0.4)


    im = axs[1].matshow(incorrect_bps, cmap='viridis', interpolation='nearest', vmin=0, vmax=1)
    axs[1].set_xticks(np.arange(0, all_diffs.shape[1], 10), np.arange(0, all_diffs.shape[1], 10))
    axs[1].set_yticks(np.arange(0, all_diffs.shape[0], 10), np.arange(0, all_diffs.shape[0], 10))
    axs[1].grid(color='gray', linestyle='--', linewidth=0.4)
    axs[1].set_title('Incorrect bps of all clusters')

    im = axs[2].matshow(missing_bps, cmap='viridis', interpolation='nearest',vmin=0, vmax=1)
    axs[2].set_xticks(np.arange(0, all_diffs.shape[1], 10), np.arange(0, all_diffs.shape[1], 10))
    axs[2].set_yticks(np.arange(0, all_diffs.shape[0], 10), np.arange(0, all_diffs.shape[0], 10))
    axs[2].grid(color='gray', linestyle='--', linewidth=0.4)
    axs[2].set_title('Missing bps of all clusters')

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    if folder:
        plt.savefig(folder + '/' + 'heatmap_diff_bp.png')
        plt.clf()
    else:
        plt.show()
        plt.clf()
    print(f'relative distance:\t{rel_distance:.2f}')
    #return result_matrix



def main():
    # simRNA_folder = '/scr/aldea/kgutenbrunner/simRNA/TBFV_sim/ALKV_xrRNA1'
    simRNA_folder = '/Users/katringutenbrunner/Desktop/MA/working/simRNA/03_14'
    # simRNA_folder = '/Users/katringutenbrunner/Desktop/MA/working/simRNA/02_27'
    #simRNA_folder = '/Users/katringutenbrunner/Desktop/MA/working/simRNA/seq2_mc_500000'
    df, [true_structure, ss_matrices] = summary_analysis(simRNA_folder) #, save=True)

    # bp_plot_vs_ts(ss_matrices, true_structure) #, simRNA_folder)
    bp_diff_plot(true_structure, ss_matrices, simRNA_folder)


if __name__ == "__main__":
    main()

