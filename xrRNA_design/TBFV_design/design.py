import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut
from tqdm import tqdm
import numpy as np

# iupac_cons =  'NNNNNNNNXXXNNNNNXGGCAGCRCRCXXNNNXXXXXXXXGYGACGGGXXXXXXXXXGGUCXXXXXXCCCGACXNNNNNXXXNNNNNNNNNNNNXXXXXXXUUYGUGAXGACCXX'
# structures = ['..(((((((..(((((((......(((((.........))))).(((((((.............))))))).)))))))..)))))))...........................',
# 		      '.....................(((...............................................................................))).........',
# 		      '.......................................................(((((((..............................................)))))))']
# #   		   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
# #    		   0        10        20        30        40        50        60        70        80        90       100       110

iupac_cons =  'NNNNNNNNXXXNNNNNXGGCAGCRCRCXXNNNXXXXXXXXGYGACGGGXXXXXXXXXGGUCXXXXXXCCCGACXNNNNNXXXNNNNNNNNNNNNXXXXXXXUUYGUGAXGACCXX'
structures = ['..(((((((..(((((((......(((((.........))))).(((((((.............))))))).)))))))..)))))))...........................',
		      '.....................(((...............................................................................))).........',
		      '.......................................................(((((((..............................................)))))))']
#   		   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
#    		   0        10        20        30        40        50        60        70        80        90       100       110

var_stem_regions = [[(27,28), (38,39)], [(48,50), (64,66)], [(55,56), (113,114)]] # start, stop -> including stop
var_loop_regions = [(9, 10), (32, 37), (51, 54), (62, 63), (79,80), (94,100)] # start, stop -> including stop

stems = {0: [[(2,9), (11, 18)], [(72, 79), (81, 88)]], 1: [(24,29), (38, 43)], 2: [(44, 51), (64, 71)]} # range
loops = {'hl2': (29,38), 'hl3': (51, 64), 'upk1': (88, 103)} # range
structure_span = {'stem': stems, 'loop': loops, 'gaps':{'beta':(27, 39), 'gamma':(48, 66)}}
# 'ss':{'beta':(8, 25), 'gamma':(26, 62)}
target_len = False # 89
target_structure = structures[0]


#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization(model_input, target_structure, start=None, steps = 1000):
    n = len(model_input.structures[0])

    # create targeted model to find weights for energy and length feature
    model_target = ir_ut.create_model_target(model_input)
    

    sampler = ir.Sampler(model_target)
    # min = 84, max = 90, mean = 86
    # need to add 2, have 2 nts in the beginning
    # 86 - 92
    sampler.set_target( 88, 4, 'totLength')
    sampler.set_target( -30, 10, 'energy')
    samples = [sampler.targeted_sample() for _ in range(1000)]
    print("Weights",{k:f'{f.weight:.3f}' for k,f in sampler.model.features.items()})

    # create model
    model = ir_ut.create_model(model_input, sampler.model.features)

    # start Monte carlo optimization
    (best_ed, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.ensemble_defect(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    start
                                    )


    # start Monte carlo optimization
    (best, best_val), sampler = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.target_frequency(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    best_ed
                                    )

    # get sequence from the MC optimization
    sample = rna.values_to_seq(best.values()[:n])
    # get structure from sample
    culled_structures = [ut.remove_positioned_gaps(sample, s) for s in structures]
    culled_seq = sample.replace('-','')

    # get target-frequency, energy of PK2, ensemble frequency, MFE and predicted secondary structure
    fc = RNA.fold_compound(culled_seq)
    fc.pf()
    (ss, mfe) = fc.mfe()
    freq = ut.target_frequency(sample, model_input)
    pk2_e = fc.eval_structure(culled_structures[2])
    ed = fc.ensemble_defect(ss)

    # print all information about best sample
    if True:
        print('\n')
        ut.margin_left('sequence:', sample, 30)
        ut.margin_left('IUPAC:', iupac_cons, 30)
        ut.margin_left('target structure:', target_structure, 30)
        print('\nculled')
        ut.margin_left('sequence:', culled_seq, 30)
        ut.margin_left('base structure:', culled_structures[0], 30)
        ut.margin_left('pk1 structure:', culled_structures[1], 30)
        ut.margin_left('pk2 structure:', culled_structures[2], 30)

        print('\nRNAFold predictions')
        ut.margin_left('target structure:', culled_structures[0], 30)
        ut.margin_left('MFE structure:', ss, 30)
        ut.margin_left('structure == MFE:', ss ==culled_structures[0], 30)
        ut.margin_left('MFE:', f'{mfe:4.2f}', 30)
        ut.margin_left('length:', len(culled_seq), 30)
        ut.margin_left('energy:', f"{RNA.energy_of_struct(culled_seq,culled_structures[0]):4.2f}", 30)
        ut.margin_left('PK2 energy:', f'{pk2_e:2.4f}', 30)
        ut.margin_left('frequency:', f'{freq:2.4f}', 30)
        ut.margin_left('ensemble defect:', f'{ed:2.4f}', 30)


def creating_samples(steps = 1000):
    model_input = ir_ut.ModelInput(structures=structures,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span,
                                target_length = target_len
                                )

    n = len(model_input.structures[0])

    # create targeted model to find weights for energy and length feature
    model_target = ir_ut.create_model_target(model_input)
    sampler = ir.Sampler(model_target)
    # min = 84, max = 90, mean = 86
    # need to add 2, have 2 nts in the beginning
    # 86 - 92
    sampler.set_target( 88, 4, 'totLength')
    sampler.set_target( -25, 15, 'energy')
    samples = [sampler.targeted_sample() for _ in range(1000)]

    # create model for sequence design
    model = ir_ut.create_model(model_input, sampler.model.features)
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(steps)]
    sequences = [rna.values_to_seq(sample.values()[:n]) for sample in samples]
    return sequences


def main():
    model_input = ir_ut.ModelInput(structures=structures,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span,
                                target_length = target_len
                                )

    mc_optimization(model_input, target_structure=target_structure, steps=100000)

if __name__ == "__main__":
    main()
