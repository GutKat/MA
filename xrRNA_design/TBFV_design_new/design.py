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
#iupac_seq =   "NNNNNNNXXXXXXGGCAGCRCRCXXXXNNNXXXXXXGYGACGGGXXXXXXXXXGGUCXXXXXXXCCCGACXXXXXXNNNNNNNNNXXXXXXXXXXXXXUUYGUXXXGAGACC"
iupac_cons = 'NNNNNNXXXXXNNNXXXGGCAGCRCNNXXNNNXXXXXXXXNNGACNNNXXXNNXXXXGGUCXXXXXXXNNNGACXXXNNNXXXXXNNNNNNNNNNNNNXXXXXXYGUXXXGACCX'
structure = ['..(((((((..(((((((......(((((.........))))).(((((((..............))))))).)))))))..)))))))..........................',
		     '.....................(((................................................................................)))........',
		     '........................................................((((((...............................................))))))']
#   		  01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#    		  0        10        20        30        40        50        60        70        80        90       100       110
var_stem_regions = [[(6,8), (82, 84)],[(14,16), (74, 76)], [(27,28), (38,39)], [(48,50), (65,67)]] # start, stop -> including stop
var_loop_regions = [(9, 10), (80,81), (32, 37), (53, 55), (62, 64), (98,103)] # start, stop -> including stop


stems = {0: [[(2,9), (11, 18)], [(73, 80), (82, 89)]], 1: [(24,29), (38, 43)], 2: [(44, 51), (65, 72)]} # range
loops = {'hl2': (29,38), 'hl3': (51, 65), 'upk1': (89, 104)} # range
structure_span = {'stem': stems, 'loop': loops}

target_len = 95
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
    ut.margin_left('target structure:', target_structure, 30)
    ut.margin_left('sequence:', sample, 30)
    if False:
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
        if False:
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
                                structure_span = structure_span,
                                target_length = target_len
                                )


    # ut.weight_testing(model_input ,target_structure, steps = 1000)
    # ut.constraint_testing(sampling_no=10000)

    model = ir_ut.create_model(model_input)
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(10)]
    possible_gaps = [i for i in range(len(iupac_cons)) if iupac_cons[i] == 'X']
    possible_gaps.sort()
    for sample in samples:
        X = sample.values()[0:len(iupac_cons)]
        Y = sample.values()[len(iupac_cons):]
        seq = rna.values_to_seq(sample.values()[0:len(iupac_cons)])
        print('X: ', X)
        print('X: ', [X[i] for i in possible_gaps])
        print('Y: ', Y)
        print(iupac_cons)
        print(seq)
        print(seq.count('-'), '\n')



if __name__ == "__main__":
    main()