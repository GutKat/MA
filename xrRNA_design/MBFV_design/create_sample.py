import infrared as ir
import infrared.rna as rna
import os
import utils as ut
import ir_utils as ir_ut
import pandas as pd
from tqdm import tqdm

# base Triples     -  +                   + -                                -     +
# db =        '..[[.((((((((((.......))))).((((((((.......{{{{{{{{.....))))))))]])))))......}}}}}}}}..'
iupac_cons =  'NNWGUCAGGCCXXXXNNNXXXXXXXXGCYACNXXXXXXXXXXXNNNXXXXXNXXXXXXXXNGUGCWGCCUGXXXXXXXXXXXNNNNN'
structures = ['.....((((((((((.......))))).((((((((....................))))))))..)))))................',
		      '..((............................................................)).....................',
		      '...........................................((((((((..........................))))))))..']
#   		   0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#    		   0        10        20        30        40        50        60        70        80

# variable regions of the design
var_stem_regions = [[(11, 14), (22, 25)], [(32, 35), (56, 59)], [(46, 50), (77, 81)]] # start, stop -> including stop
var_loop_regions = [(18, 21), (36, 42), (52,55), (71, 76)] # start, stop -> including stop

# alpha (0), beta (1), gamma (2) structure domains
stems = {0: [(5, 10), (66, 71)], 1: [(10, 15), (22, 27)], 2:[(28, 36), (56, 64)]} # range
loops = {'hl2': (15,22), 'hl3': (36, 56), 'upk1': (71, 77)} # range

# all the different regions of the structure saved
# ss and gaps are used for limiting the length of beta and gamma
# ss is the whole range of the beta and gamma structure, gaps is the position of the first and last possible gap within the structure
structure_span = {'stem': stems, 'loop': loops, 'ss':{'beta':(11, 27), 'gamma':(28, 64)}, 'gaps':{'beta':(12, 25), 'gamma':(32, 59)}}

# target_len can be an int, a range, or False
# int = excat sequence length, range = sequence length is within range, False = no length constraint target_len = False #False # [52, 70]
target_len = False # [53, 78]
target_structure = structures[0]


#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization_testing(model_input, steps = 100000):

    n = len(model_input.structures[0])
    model_target = ir_ut.create_model_target(model_input)
    sampler = ir.Sampler(model_target)
    # biological data ranges from 53 - 63 (mean=60)
    # need to add 2 because we have 2 unpaired nt in beginning --> 55 - 65, mean =62
    sampler.set_target(59, 7, 'totLength')
    sampler.set_target( -18, 12, 'energy')
    # sampler.set_target( -4, 3, 'energy_pk2')
    samples = [sampler.targeted_sample() for _ in range(5000)]
    # create model for sequence design
    model = ir_ut.create_model(model_input, sampler.model.features)
    sampler = ir.Sampler(model)
    start = sampler.sample()
    before_optimization = rna.values_to_seq(start.values()[:n])

    (best_ed, best_val) = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.ensemble_defect(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    start
                                    )
    # start Monte carlo optimization
    (best, best_val) = ut.mc_optimize(model,
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
    folder_name = '/scr/aldea/kgutenbrunner/working/xrRNA_design/MBFV_design/seqs/before_after_opt/'
    csv_file = 'designs_before_after_opt'
    model_input = ir_ut.ModelInput(structures=structures,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span,
                                target_length = target_len
                                )
    sequences_before = []
    sequences_after = []
    num_tasks = 1000

    # Automatically determine the number of available CPU cores to maximize parallelism
    print(multiprocessing.cpu_count())
    max_workers =multiprocessing.cpu_count() #8  # or set manually, e.g., 4, 8, etc.
    # Using ProcessPoolExecutor to run as many tasks in parallel as possible
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submitting the mc_optimization_testing function for parallel execution
        futures = [executor.submit(mc_optimization_testing, model_input, 100000) for _ in range(num_tasks)]

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


