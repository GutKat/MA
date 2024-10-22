import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut
import time

# base Triples     -  +                   + -                                -     + 
db =          '..[[.((((((((((.......))))).((((((((.......{{{{{{{{.....))))))))]])))))......}}}}}}}}..'
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
target_structure = structures[0]

# monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_optimization(model_input, target_structure, start=None, steps=100000):
    # create model
    model = ir_ut.create_model(model_input)
    n = len(model_input.structures[0])

    # start Monte carlo optimization
    tic = time.time()
    sampler = ir.Sampler(model)
    # biological data ranges from 50 - 62 (mean=57)
    # need to add 2 because we have 2 unpaired nt in beginning --> 52 - 64, mean = 59
    sampler.set_target(59, 7, 'totLength')
    sampler.set_target( -18, 12, 'energy', -1)
    (best, bestval),neighbors_top10, sampler = ut.neighbor_optimize_top10(sampler,
                                                              model,
                                                              model_input,
                                                              lambda sequence: ut.ensemble_defect(sequence, model_input),
                                                              round(steps),
                                                              0.1,
                                                              0.15,
                                                              0.05,
                                                              rna.HammingDist,
                                                              start,
                                                              )

    (best, bestval),neighbors_top10, sampler = ut.neighbor_optimize_top10(sampler,
                                                              model,
                                                              model_input,
                                                              lambda sequence: ut.target_frequency(sequence, model_input),
                                                              round(steps),
                                                              0.1,
                                                              0.15,
                                                              0.05,
                                                              rna.HammingDist,
                                                              neighbors_top10,
                                                              )

    tac = time.time()
    # get sequence from the MC optimization
    sample = rna.values_to_seq(best.values()[:n])
    # get structure from sample
    culled_structures = [ut.remove_positioned_gaps(sample, s) for s in structures]
    culled_seq = sample.replace('-', '')

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
        ut.margin_left('target structure:', db, 30)
        print('\nculled')
        ut.margin_left('sequence:', culled_seq, 30)
        ut.margin_left('base structure:', culled_structures[0], 30)
        ut.margin_left('pk1 structure:', culled_structures[1], 30)
        ut.margin_left('pk2 structure:', culled_structures[2], 30)

        print('\nRNAFold predictions')
        ut.margin_left('target structure:', culled_structures[0], 30)
        ut.margin_left('MFE structure:', ss, 30)
        ut.margin_left('structure == MFE:', ss == culled_structures[0], 30)
        ut.margin_left('MFE:', f'{mfe:4.2f}', 30)
        ut.margin_left('length:', len(culled_seq), 30)
        ut.margin_left('energy:', f"{RNA.energy_of_struct(culled_seq, culled_structures[0]):4.2f}", 30)
        ut.margin_left('PK2 energy:', f'{pk2_e:2.4f}', 30)
        ut.margin_left('frequency:', f'{freq:2.4f}', 30)
        ut.margin_left('ensemble defect:', f'{ed:2.4f}', 30)


    print(f'time: {round((tac-tic)/60,3)}')

def main():
    model_input = ir_ut.ModelInput(structures=structures,
                                   anti_structures=[],
                                   iupac=iupac_cons,
                                   var_stem_regions=var_stem_regions,
                                   var_loop_regions=var_loop_regions,
                                   structure_span=structure_span,
                                   )

    mc_optimization(model_input, target_structure=structures, steps=1)


if __name__ == "__main__":
    main()
