import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut

# base Triples   -  +                   + -                                -     + 
# db =        '[[.((((((((((.......))))).((((((((.......{{{{{{{{.....))))))))]])))))......}}}}}}}}..'
iupac_cons =  'WGUCAGGCCXXXXNNNXXXXXXXXGCYACNXXXXXXXXXXXNNNXXXXXNXXXXXXXXNGUGCWGCCUGXXXXXXXXXXXNNNNN'
structures = ['...((((((((((.......))))).((((((((....................))))))))..)))))................',
		      '((............................................................)).....................',
		      '.........................................((((((((..........................))))))))..']
#   		   0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#    		   0        10        20        30        40        50        60        70   

# variable regions of the design
var_stem_regions = [[(9, 12), (20, 23)], [(30, 33), (54, 57)], [(44, 48), (75, 79)]] # start, stop -> including stop
var_loop_regions = [(16, 19), (34, 40), (50,53),  (69, 74)] # start, stop -> including stop

# alpha (0), beta (1), gamma (2) structure domains
stems = {0: [(3, 8), (64, 69)], 1: [(8, 13), (20, 25)], 2:[(26, 34), (54, 62)]} # range
loops = {'hl2': (13,20), 'hl3': (34, 54), 'upk1': (69, 75)} # range

# all the different regions of the structure saved
# ss and gaps are used for limiting the length of beta and gamma
# ss is the whole range of the beta and gamm structure, gaps is the position of the first and last possible gap within the structure
structure_span = {'stem': stems, 'loop': loops, 'ss':{'beta':(8, 25), 'gamma':(26, 62)}, 'gaps':{'beta':(10, 23), 'gamma':(30, 57)}}

# beta_min = 7; beta_max = 15
# gamma_min = 20; gamma_max = 23

# target_len can be an int, a range, or False
# int = excat sequence length, range = sequence length is within range, False = no length constraint target_len = False #False # [52, 70]
target_len = [52, 70] #[52, 70] False
target_structure = structures[0]


#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization(model_input, target_structure, start=None, steps = 100000):
    # create model
    model = ir_ut.create_model(model_input)
    n = len(model_input.structures[0])

    # start Monte carlo optimization
    (best, best_val), sampler = result = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence: ut.objective_function(sequence, model_input),
                                    steps,
                                    0.01,
                                    start
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
        ut.margin_left('folding:', ss, 30)
        ut.margin_left('MFE:', f'{mfe:4.2f}', 30)
        ut.margin_left('length:', len(culled_seq), 30)
        ut.margin_left('energy:', f"{RNA.energy_of_struct(culled_seq,culled_structures[0]):4.2f}", 30)
        ut.margin_left('PK2 energy:', f'{pk2_e:2.4f}', 30)
        ut.margin_left('frequency:', f'{freq:2.4f}', 30)
        ut.margin_left('ensemble defect:', f'{ed:2.4f}', 30)
        ut.margin_left('objective funtion:', f'{best_val:2.4f}', 30)




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


######################################
# for testing the objective function #
######################################

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
