from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
from collections import namedtuple
import random
import math
import ir_utils as ir_ut
import itertools as it



def target_frequency(sequence, target):
    ss = remove_positioned_gaps(sequence, target)
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)


def count_gaps(sequence, region):
    return sequence[region[0]:region[1]].count('-')


def optimization_function(sequence, model_input):
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']

    st1 = count_gaps(sequence, stems[0][0]) + count_gaps(sequence, stems[0][1])
    st2 = count_gaps(sequence, stems[1][0]) + count_gaps(sequence, stems[1][1])
    st3 = count_gaps(sequence, stems[2][0]) + count_gaps(sequence, stems[2][1])#
    hl2 = count_gaps(sequence, loops['hl2'])
    hl3 = count_gaps(sequence, loops['hl3'])
    upk1 = count_gaps(sequence, loops['upk1'])
    print(st1, st2, st3, hl2, hl3, upk1)

    target = model_input.structures[0]
    fc_pr = target_frequency(sequence, target)
    print(fc_pr)
    return fc_pr


def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)

    


def mc_optimize(model, model_input, objective, steps, temp, start=None):
    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start

    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1/len(cc) for cc in ccs]
    for i in (range(steps)): #tqdm
        cc = random.choices(ccs,weights)[0]
        new = sampler.resample(cc, cur)
        new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
        newval = objective(new_seq)
        if (newval >= curval or random.random() <= math.exp((newval-curval)/temp)):
            cur, curval = new, newval
            if curval > bestval:
                #print(curval)
                best, bestval = cur, curval
        print(f'{i} step')
    return (best, bestval)


def weight_testing(model_input, target_structure, steps = 1000):
    model = ir_ut.create_model(model_input)

    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']

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
        sample = rna.values_to_seq(sample.values()[:len(model_input.structures[0])])
        print(sample)
        print(model_input.structures[0])



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


def constraint_testing(sampling_no = 100):
    #              012345678901
    test_ss =    ['((((..(())..))))']
    iupac_cons =  'NNXXNNXXXXNN'
    var_stem_regions = [[(2, 3), (8, 9)], [(6,7), (8,9)]]  # start, stop -> including stop
    var_loop_regions = [(4, 5)] # start, stop -> including stop
    structure_span = {}
    model_input = ir_ut.ModelInput(structures=test_ss, #structure_w_gaps
                                anti_structures=[],
                                iupac=iupac_cons,
                                var_stem_regions=var_stem_regions,
                                var_loop_regions=var_loop_regions,
                                structure_span = structure_span
                                )
    model = ir_ut.create_model(model_input)
    sampler = ir.Sampler(model)
    samples = [sampler.sample() for _ in range(sampling_no)]
    tests_failed = False
    failed_sample = None
    for sample in samples:
        sample = rna.values_to_seq(sample.values()[:len(model_input.structures[0])])
        # check if sample passes constraint test
        pass_test = ir_ut.constraint_testing([iupac_cons, var_loop_regions, var_stem_regions])
        if not pass_test:
            print('!' * 50)
            print(f'Test not passed!')
            print('!' * 50)
            tests_failed = True
            failed_sample = sample
        print(sample)
        print(test_ss[0])
    if tests_failed:
        print('!' * 50)
        print(f'Test not passed!')
        print(failed_sample)
        print('!' * 50)
    else:
        print('test passed - all samples fulfill constraints')


def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")