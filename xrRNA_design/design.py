import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut

'''
structure overview

S-I = ST-I + IL-I = (22 + 0) - (28 + 4) = 22 - 32
--> could make it variable where the internal loop is since there seems to be no consensus -> but for now, keep it in the middle
S-II = ST-II + HL-II = (6 + 3) - (10 + 9) = 9 - 19
S-III = ST-III + HL-III = (8 + 4) - (14 + 12) = 12 - 26
--> longer because I added a variable PK2 -> HL-III could be now 14 nt long, which is not in the data, but it was easier this way
PKs = uPK + PK1 + PK1bPK2 + PK2 = (8 + 3 + 0 + 4) - (15 + 3 + 2 + 6) = 15 - 26
'''

# cons =
iupac_cons = 'NNNNXXXXXNNNXXXGGCAGCRCNNXXNNNXXXXXXXXNNGACNNNXXXNNXXXXGGUCXXXXXXXNNNGACXXXNNNXXXXXNNNNNNNNNNNNNXXXXXXUGYXXXGACCX'
structure = ['(((((((..(((((((......(((((.........))))).(((((((..............))))))).)))))))..)))))))..........................',
		     '...................(((................................................................................)))........',
		     '......................................................((((((...............................................))))))']
#   		  01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#    		  0        10        20        30        40        50        60        70        80        90       100       110
len('(((((((..(((((((......(((((.........))))).(((((((..............))))))).)))))))..)))))))..........................')
var_stem_regions = [[(4,6), (80, 82)],[(12,14), (72, 74)], [(25,26), (36,37)], [(46,48), (63,65)]] # start, stop -> including stop
var_loop_regions = [(7, 8), (78,79), (30, 35), (51, 53), (60, 62), (96,101)] # start, stop -> including stop

stems = {0: [[(0,7), (9, 16)], [(71, 78), (80, 87)]], 1: [(22,27), (36, 41)], 2: [(42, 49), (63, 70)]} # range
loops = {'hl2': (27,36), 'hl3': (49, 63), 'upk1': (87, 102)} # range
structure_span = {'stem': stems, 'loop': loops}

target_len = 87
target_gc =  0.58
target_energy = -33
target_structure = structure[0]

#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_negative_optimization(model_input, target_structure, start=None, steps = 100000):
    model = ir_ut.create_model(model_input)
    n = len(model_input.structures[0])

    (best, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.optimization_function(sequence, model_input), # add relations of strucutre into functions
                                    steps,
                                    0.01,
                                    start
                                   )


    sample = rna.values_to_seq(best.values()[:n])
    samples = [sampler.sample() for _ in range(10)]

    culled_structure = ut.remove_positioned_gaps(sample, target_structure)
    culled_seq = sample.replace('-','')

    fc = RNA.fold_compound(culled_seq)
    fc.pf()
    (ss, mfe) = fc.mfe()
    freq = ut.target_frequency(sample, target_structure)

    if True:
        print('\n')
        ut.margin_left('target structure:', target_structure, 30)
        ut.margin_left('sequence:', sample, 30)
        ut.margin_left('cullled structure:', culled_structure, 30)
        ut.margin_left('cullled sequence:', culled_seq, 30)
        print('\nRNAFold predictions')
        ut.margin_left('folding:', ss, 30)
        ut.margin_left('MFE:', f'{mfe:4.2f}', 30)
        ut.margin_left('length:', len(culled_seq), 30)
        ut.margin_left('energy:', f"{RNA.energy_of_struct(culled_seq,culled_structure):4.2f}", 30)
        ut.margin_left('frequency:', f'{freq:2.4f}', 30)
        ut.margin_left('objective funtion:', f'{best_val:2.4f}', 30)
        if True:
            with open("/scr/aldea/kgutenbrunner/working/xrRNA_design/TBFV_design/seqs/seq.out", "w") as file:
                file.write(sample)
    # print('\n')
    # for sample in samples:
    #     sample = rna.values_to_seq(sample.values()[:len(model_input.structures[0])])
    #     culled_structure = ut.remove_positioned_gaps(sample, target_structure)
    #     culled_seq = sample.replace('-','')

    #     fc = RNA.fold_compound(culled_seq)
    #     fc.pf()
    #     (ss, mfe) = fc.mfe()
    #     freq = ut.target_frequency(sample, target_structure)
    #     # ut.margin_left('frequency:', f'{freq:2.4f}', 15)




def main():

    model_input = ir_ut.ModelInput(structures=structure, #structure_w_gaps
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span
                                )


    # ut.weight_testing(model_input ,target_structure, steps = 1000)
    # ut.constraint_testing(sampling_no=10000)
    
    mc_negative_optimization(model_input, target_structure= target_structure, steps=100000)


if __name__ == "__main__":
    main()