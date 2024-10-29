import infrared as ir
import infrared.rna as rna
import os
import utils as ut
import ir_utils as ir_ut
import pandas as pd
from tqdm import tqdm

iupac_cons =  'NNNNNNNNXXXNNNNXXGGCAGCRCRCXXNNNXXXXXXXXGYGACGGGXXXXXXXXGGUCXXXXXXCCCGACXXNNNNXXXNNNNNNNNNNNNXXXXXXXUUYGUGAXGACCXX'
structures = ['..(((((((..(((((((......(((((.........))))).(((((((............))))))).)))))))..)))))))...........................',
		      '.....................(((..............................................................................))).........',
		      '......................................................(((((((..............................................)))))))']
#   		   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
#    		   0        10        20        30        40        50        60        70        80        90       100       110

var_stem_regions = [[(15,16), (72,73)], [(27,28), (38,39)], [(48,50), (63,65)], [(54,55), (112,113)]] # start, stop -> including stop
var_loop_regions = [(9, 10), (32, 37), (51, 53), (61, 62), (78,79), (93,99)] # start, stop -> including stop

gaps = {'alpha': (8, 16, 72, 80), 'beta':(27, 39), 'gamma':(48, 65), 'gamma_hl': (51, 62), 'gamma_stem': (8, 50, 63, 65)}
# 'ss':{'beta':(8, 25), 'gamma':(26, 62)}
target_len = False # 89
target_structure = structures[0]


#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization_testing(model_input, steps = 100000):

    n = len(model_input.structures[0])

    # create model
    model = ir_ut.create_model(model_input)
    sampler = ir.Sampler(model)
    start = sampler.sample()
    before_optimization = rna.values_to_seq(start.values()[:n])

    (best_ed, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.ensemble_defect(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    start
                                    )
    print(best_ed)
    # start Monte carlo optimization
    (best, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.target_frequency(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    best_ed
                                    )

    # get sequence from the MC optimization
    after_optimization = rna.values_to_seq(best.values()[:n])
    return before_optimization, after_optimization


def save_csv_with_suffix(df, folder_path, file_name):
    # Ensure the folder path exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Construct the full path of the initial CSV file
    full_path = os.path.join(folder_path, f"{file_name}.csv")

    # Initialize a counter for the file suffix
    counter = 1

    # Check if the file already exists
    while os.path.isfile(full_path):
        # Modify the file name by appending _1, _2, etc.
        full_path = os.path.join(folder_path, f"{file_name}_{counter}.csv")
        counter += 1

    # Save the DataFrame to the determined path
    df.to_csv(full_path, index=False)
    print(f"File saved as: {full_path}")


import concurrent.futures
import multiprocessing

if __name__ == "__main__":
    folder_name = '/scr/aldea/kgutenbrunner/working/xrRNA_design/TBFV_design/data/seqs/before_after_opt/29_10_24/'
    csv_file = 'designs_before_after_opt'
    model_input = ir_ut.ModelInput(structures=structures,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                gaps = gaps,
                                target_length = target_len
                                )
    sequences_before = []
    sequences_after = []
    num_tasks = 10


    # Automatically determine the number of available CPU cores to maximize parallelism
    max_workers =multiprocessing.cpu_count() #8  # or set manually, e.g., 4, 8, etc.
    # Using ProcessPoolExecutor to run as many tasks in parallel as possible
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submitting the mc_optimization_testing function for parallel execution
        futures = [executor.submit(mc_optimization_testing, model_input, 100) for _ in range(num_tasks)]

        # Collecting the results as they are completed
        for future in tqdm(concurrent.futures.as_completed(futures), total=num_tasks):
            before_opt, after_opt = future.result()  # Get result from the future
            sequences_before.append(before_opt)
            sequences_after.append(after_opt)

    # Combining before and after sequences into a DataFrame
    sequences = {'before': sequences_before, 'after': sequences_after}
    df = pd.DataFrame(sequences)
    print(df)

    # Saving the DataFrame to a CSV file
    save_csv_with_suffix(df, folder_name, csv_file)


