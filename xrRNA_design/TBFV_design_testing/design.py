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
iupac_cons = 'NNNXXXNXX'
structure = ['.........']
#   		  01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012
#    		  0        10        20        30        40        50        60        70        80        90       100       110

loops = {'hl2': (0,3),  'upk1': (7, 9)} # range
structure_span = { 'loop': loops}

target_len = 7
target_gc =  0.58
target_energy = -33
target_structure = structure[0]

#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_negative_optimization(model_input, target_structure, start=None, steps = 100000):
    model = ir_ut.create_model(model_input)
    n = len(model_input.structures[0])

    (best, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.target_frequency(sequence, model_input), # add relations of strucutre into functions
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
    freq = ut.target_frequency(sample, model_input)
    ut.margin_left('target structure:', target_structure, 30)
    ut.margin_left('sequence:', sample, 30)

def main():
    model_input = ir_ut.ModelInput(structures=structure, #structure_w_gaps
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=[],
                                var_loop_regions=[],
                                structure_span = structure_span,
                                target_length = target_len
                                )


    model = ir_ut.create_model(model_input)
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(1)]
    possible_gaps = [i for i in range(len(iupac_cons)) if iupac_cons[i] == 'X']
    n_y = len(possible_gaps)
    possible_gaps.sort()
    for sample in samples:
        model.eval_fat
        print(sample.values())
        X = sample.values()[0:len(iupac_cons)]
        Y = sample.values()[len(iupac_cons):len(iupac_cons)+n_y]
        Z = sample.values()[len(iupac_cons)+n_y:]
        seq = rna.values_to_seq(sample.values()[0:len(iupac_cons)])
        print('X: ', X)
        print('X: ', [X[i] for i in possible_gaps])
        print('Y: ', Y)
        print('Z: ', Z)
        print(iupac_cons)
        print(seq)
        print('gaps: ',seq.count('-'), '\n')


if __name__ == "__main__":
    main()