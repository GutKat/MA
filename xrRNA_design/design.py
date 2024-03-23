import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
import ir_utils as ir_ut


#       *(((((((((...GCR(((..))).(((....GGUC..))).))))))))).........UGYGACC
#       *(((((((((...[[[(((..))).(((....{{{{..))).))))))))).........]]]}}}}         minimal structure
# ST        9-14           3-5          3-5          9-14
# IL/HL     0-2            2-4        10-12         0-2       9-16    0/2
# PKs                 3                  4-7                        3  4-7

# maximal structure:
#       *((((((..((((((((...[[[(((((....))))).(((((....{{{{{{{....))))).))))))))..))))))................]]]..}}}}}}}         
#       '((((XXXX((((XXXGGCAGCRC((XX..XXXX))GAC((XX....XXXGGUC..XXXX))GACXXX))))XXXX))))..........XXXXXXUGYXXGACCXXX''
#        012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#        0        10         20       30        40        50        60        70        80        90       100       110


#                           ST-I A               ST-I B           ST - II          ST-III                  PK2
var_stem_regions = [[(4,5), (73, 74)],[(12,14), (64, 66)], [(25,26), (31,32)], [(40,41), (57,58)], [(46,48), (104,106)]]
var_loop_regions = [(6,7), (29,30), (55,56), (71,72), (89,94), (98,99)]

stems = {0: [(0,16), (63, 79)], 1: [(22,27), (31, 36)], 2: [(37, 42), (57, 62)]}
loops = {'hl2': (27,31), 'hl3': (42, 57), 'upk1': (79, 94)}


ss ='((((((..((((((((...[[[(((((....))))).(((((....{{{{{{{....))))).))))))))..))))))................]]]..}}}}}}}'

# cons =            '((((XXXX((((XXXGGCAGCRC((XX..XXXX))GAC((XX....XXXGGUC..XXXX))GACXXX))))XXXX))))..........XXXXXXUGYXXGACCXXX'
iupac_cons =        'NNNNXXXXNNNNXXXGGCAGCRCNNXXNNXXXXNNGACNNXXNNNNXXXGGUCNNXXXXNNGACXXXNNNNXXXXNNNNNNNNNNNNNNXXXXXXUGYXXGACCXXX'
structure =        ['((((((..((((((((......(((((....))))).(((((...............))))).))))))))..))))))............................',
                    '...................(((.........................................................................))).........',
                    '.................................................((((...............................................))))...']
structure_w_gaps = ['((((XX..((((XXX(......(((XX....XX))).(((XX...............XX))).)XXX))))..XX))))............................',
                    '...................(((.........................................................................))).........',
                    '..............................................XXX((((...............................................))))XXX']
#                    012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#                    0        10         20       30        40        50        60        70        80        90       100       110


target_len = 94
target_gc =  0.58
target_structure = structure[0]

#monte carlo optimization of the sequence design - objective function is frequency of target structure
def mc_negative_optimization(model_input, target_structure, start=None, steps = 100000):
    model = ir_ut.create_model(model_input)
    

    n = len(model_input.structures[0])

    best, best_val = ut.mc_optimize(model,
                                    model_input,
                                    lambda sequence,target_structure=target_structure: ut.target_frequency(sequence, target_structure),
                                    steps,
                                    0.01,
                                    start
                                   )
    
    sample = rna.values_to_seq(best.values()[:n])

    culled_structure = ut.remove_positioned_gaps(sample, target_structure)
    culled_seq = sample.replace('-','')

    fc = RNA.fold_compound(culled_seq)
    fc.pf()
    (ss, mfe) = fc.mfe()
    print(culled_structure)
    print(f"{sample},\n{ss}:\n{mfe:4.2f}, {RNA.energy_of_struct(culled_seq,culled_structure):4.2f}, {best_val:2.4f}")


### Testing model creation and structure lengths ###
def testing(structures,target_structure, var_loop_regions, var_stem_regions, iupac = '', steps = 1000):


    model_input = ir_ut.ModelInput(structures=structures,
                                    anti_structures=[],
                                    iupac=iupac,
                                    var_stem_regions=var_stem_regions,
                                    var_loop_regions=var_loop_regions,
                          )
    model = ir_ut.create_model(model_input)

    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(steps)]

    st1_lens = []
    hl2_lens = []
    st2_lens = []
    hl3_lens = []
    st3_lens = []
    upk1_lens = []
    tot_lens = []
    for sample in samples:
        sample = rna.values_to_seq(sample.values()[:len(structures[0])])
        print(sample)
        print(structure[0])



        st1 = len(sample[slice(*stems[0][0])]+sample[slice(*stems[0][1])].replace('-', ''))

        hl2 = len(sample[slice(*loops['hl2'])].replace('-', ''))
        st2 = len(sample[slice(*stems[1][0])]+sample[slice(*stems[1][1])].replace('-', ''))


        hl3 = len(sample[slice(*loops['hl3'])].replace('-', ''))
        st3 = len(sample[slice(*stems[2][0])]+sample[slice(*stems[2][1])].replace('-', ''))

        upk1 = len(sample[slice(*loops['upk1'])].replace('-', ''))


        st1_lens.append(st1)
        
        hl2_lens.append(hl2)
        st2_lens.append(st2)

        hl3_lens.append(hl3)
        st3_lens.append(st3)
        upk1_lens.append(upk1)

        tot_lens.append(len(sample.replace('-', '')))

    print('-' * 150)

    print('ST1\t', sum(st1_lens) / len(st1_lens))

    print('HL2\t', sum(hl2_lens) / len(hl2_lens))
    print('ST2\t', sum(st2_lens) / len(st2_lens))

    print('HL3\t', sum(hl3_lens) / len(hl3_lens))
    print('ST3\t', sum(st3_lens) / len(st3_lens))

    print('uPK1\t', sum(upk1_lens) / len(upk1_lens))

    print('totLen\t', sum(tot_lens) / len(tot_lens))



def main():

    model_input = ir_ut.ModelInput(structures=structure_w_gaps,
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                )


    testing(structure_w_gaps,target_structure, var_loop_regions, var_stem_regions, iupac=iupac_cons)
    #mc_negative_optimization(model_input, target_structure= target_structure, steps=1000)


if __name__ == "__main__":
    main()