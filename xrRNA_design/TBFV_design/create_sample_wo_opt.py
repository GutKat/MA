import design
import RNA
import pandas as pd
import os
import matplotlib.pyplot as plt

# true if working from tbi, false is working from home
tbi = True

if tbi:
    base_folder = '/scr/aldea/kgutenbrunner/working'
else:
    base_folder = '/Users/katringutenbrunner/Desktop/MA/working'

def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def target_frequency(sequence, structures):
    ss = remove_positioned_gaps(sequence, structures[0])
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)



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


csv_file = os.path.join(base_folder, 'xrRNA_design/TBFV_design/data/seqs/multiple_designs_wo_opt.csv')
sequences = design.creating_samples()
df = pd.DataFrame(sequences, columns =['sequence'])
print(df)
df.to_csv(csv_file)
