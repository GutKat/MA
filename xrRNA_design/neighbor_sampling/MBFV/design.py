import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut
import numpy as np

# base Triples     -  +                   + -                                -     + 
# db =        '..[[.((((((((((.......))))).((((((((.......{{{{{{{{.....))))))))]])))))......}}}}}}}}..'
iupac_cons =  'NNWGUCAGGCCXXXXNNNXXXXXXXXGCYACNXXXXXXXXXXXNNNXXXXXNXXXXXXXXNGUGCWGCCUGXXXXXXXXXXXNNNNN'
structures = ['.....((((((((((.......))))).((((((((....................))))))))..)))))................',
		      '..((............................................................)).....................',
		      '...........................................((((((((..........................))))))))..']
#   		   0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#    		   0        10        20        30        40        50        60        70        80

var_stem_regions = [[(11, 14), (22, 25)], [(32, 35), (56, 59)], [(46, 50), (77, 81)]] # start, stop -> including stop
var_loop_regions = [(18, 21), (36, 42), (52,55), (71, 76)] # start, stop -> including stop
stems = {0: [(5, 10), (66, 71)], 1: [(10, 15), (22, 27)], 2:[(28, 36), (56, 64)]} # range
loops = {'hl2': (15,22), 'hl3': (36, 56), 'upk1': (71, 77)} # range
structure_span = {'stem': stems, 'loop': loops, 'ss':{'beta':(10, 27), 'gamma':(28, 64)}, 'gaps':{'beta':(12, 25), 'gamma':(32, 59)}}
target_structure = structures[0]


#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization(model_input, target_structure, start=None, steps = 100000):
    n = len(model_input.structures[0])
    # create targeted model to find weights for energy and length feature
    model_target = ir_ut.create_model_target(model_input)
    
    sampler = ir.Sampler(model_target)
    sampler.set_target( 62, 10, 'totLength')
    sampler.set_target( -18, 12, 'energy')
    # sampler.set_target( -4, 3, 'energy_pk2')
    sampler.targeted_sample() for _ in range(5000)
    print("Weights",{k:f'{f.weight:.3f}' for k,f in sampler.model.features.items()})

    # create model for sequence design
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

    print("Weights",{k:f'{f.weight:.3f}' for k,f in sampler.model.features.items()})
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


####################################################
# for testing the objective function and diversity #
####################################################


def testing(steps):
    model_input = ir_ut.ModelInput(structures=structures,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span,
                                target_length = target_len
                                )
    # create model
    model = ir_ut.create_model(model_input)
    n = len(model_input.structures[0])

    # start Monte carlo optimization
    (best_ed, best_val), sampler = result = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.ensemble_defect(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    None
                                    )


        # start Monte carlo optimization
    (best, best_val), sampler = result = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.target_frequency(sequence, model_input),
                                    round(steps * 0.5),
                                    0.01,
                                    best_ed
                                    )

    return rna.values_to_seq(best.values()[:n])


# #monte carlo optimization of the sequence design - objective function is frequency of target structure
# # used for testing the objective function
# def mc_optimization_testing(model_input, target_structure, start=None, steps = 100000, obj=1):
#     model = ir_ut.create_model(model_input)
#     n = len(model_input.structures[0])
#     if obj==0:
#         result = ut.mc_optimize(model,
#                                         model_input,
#                                         lambda sequence: ut.objective_function0(sequence, model_input),
#                                         steps,
#                                         0.01,
#                                         start
#                                         )
#     elif obj==1:
#         result = ut.mc_optimize(model,
#                                         model_input,
#                                         lambda sequence: ut.objective_function1(sequence, model_input),
#                                         steps,
#                                         0.01,
#                                         start
#                                         )
#     elif obj==2:
#         result = ut.mc_optimize(model,
#                                         model_input,
#                                         lambda sequence: ut.objective_function2(sequence, model_input),
#                                         steps,
#                                         0.01,
#                  start
#                                         )

#     return rna.values_to_seq(best.values()[:n])



# def testing(steps, objective_function):
#     model_input = ir_ut.ModelInput(structures=structures,
#                                 anti_structures=[],
#                                 iupac=iupac_cons,
#                                 var_stem_regions=var_stem_regions,
#                                 var_loop_regions=var_loop_regions,
#                                 structure_span = structure_span,
#                                 target_length = target_len
#                                 )

#     sample = mc_optimization_testing(model_input, target_structure=target_structure, steps=steps, obj=objective_function)
#     culled_structures = [ut.remove_positioned_gaps(sample, s) for s in structures]
#     culled_seq = sample.replace('-','')
#     fc = RNA.fold_compound(culled_seq)
#     fc.pf()
#     (ss, mfe) = fc.mfe()
#     freq = fc.pr_structure(culled_structures[0])
#     pk2_e = fc.eval_structure(culled_structures[2])
#     ed = fc.ensemble_defect(ss)
#     return sample, freq, ed, pk2_e
