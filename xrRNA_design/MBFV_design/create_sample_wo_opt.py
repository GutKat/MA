# %%
import design
import pandas as pd
from tqdm import tqdm
import RNA
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


def calculate_bpps(sequences):
    bps = []

    for seq in sequences:
        fc = RNA.fold_compound(seq)
        fc.pf()
        # (ss, mfe) = fc.mfe()
        BPP = fc.bpp()
        BPP = np.triu(BPP) + np.triu(BPP, 1).T
        bps.append(np.array(BPP))

    bps = np.sum(bps, axis=0)
    bps = bps[1:, 1:]
    bps /= len(sequences)
    return bps


def calculate_bpps_region(bpps, region):
    return np.array([bpps[row, col] for row, col in zip(range(*region[0]), range(*region[1], -1))])


def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def append_to_csv(csv_file, new_df):
    # Read the existing CSV into a DataFrame
    try:
        existing_df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"File not found. A new file will be created at {csv_file}")
        existing_df = pd.DataFrame()

    # Append the new data to the existing DataFrame
    updated_df = pd.concat([existing_df, new_df], ignore_index=True)

    # Save the updated DataFrame back to the CSV
    updated_df.to_csv(csv_file, index=False)

    return updated_df
# %%



def main():
    csv_file = '/scr/aldea/kgutenbrunner/github/MA/thesis/data/MBFV_multiple_designs_wo_opt_new.csv'
    sequences = design.creating_samples(1000)
    df = pd.DataFrame(sequences, columns=['sequence'])
    #print(df)
    #df.to_csv(csv_file)

    # calculate the mfe structure for each sequence and create a heatmap of it

    regions = {
        'alpha': [(5, 10), (70, 65)],
        'beta': [(10, 15), (26, 21)],
        'gamma': [(28, 36), (63, 55)],
        'PK1': [(2, 4), (65, 63)],
        'PK2': [(43, 57), (84, 76)],
    }

    bps = calculate_bpps(sequences)

    bpps_sum_gamma = sum(calculate_bpps_region(bps, regions['gamma'] ))
    bpps_sum_PK2 = sum(calculate_bpps_region(bps, regions['PK2'] ))

    bpps_range = [4/6, 6/4]
    print(bpps_sum_gamma,bpps_sum_PK2)
    print(bpps_sum_gamma/bpps_sum_PK2)

    # get base pairing heatmap of a matrix of base pair probabilities
    fig, ax = plt.subplots()
    cimg = ax.imshow(bps, cmap='viridis', interpolation='nearest')
    ax.set_xticks(np.arange(0, bps.shape[1], 5), np.arange(0, bps.shape[1], 5))
    ax.set_yticks(np.arange(0, bps.shape[0], 5), np.arange(0, bps.shape[0], 5))
    ax.grid(color='gray', linestyle='--', linewidth=0.4)
    cbar = plt.colorbar(cimg, ax=ax)

    cbar_ax = cbar.ax  # Access the color bar's axes
    max_value = np.max(bps)
    cbar_ax.axhline(y=max_value, color='red', linestyle='--')  # Add a horizontal line
    plt.show()


   high_value = np.where(bps > 0.1)
   print("High-value:")
   for x,y in zip(high_value[0], high_value[1]):
       print(x,y)






if __name__ == "__main__":
    main()
    #sequences = creating_samples(5)




# %%
