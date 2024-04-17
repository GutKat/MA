from tqdm import tqdm
import infrared as ir
import RNA
import infrared.rna as rna
from collections import namedtuple
import random
import math
import ir_utils as ir_ut
import matplotlib.pyplot as plt
import numpy as np


def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")


target_len = 94
target_gc =  0.58
target_energy = -33
relations = {
        'st2_hl2': [0.81, 0.3, 1.75],
        'st3_hl3': [0.48, 0.3, 1.75],
        'sI_sII': [1.56, 0.66, 2.29],
        'sII_sIII': [0.63, 0.4, 1],
        'sI_sIII': [0.91, 0.46, 1.2]
    }
length_fraction = 0.0
relation_fraction = 0.05 
plotting_steps = 100


def number_in_range(x, range_):
    if range_[0] <= x <= range_[1]:
        return True
    return False


def error(x, y):
    return abs(x-y)


def remove_positioned_gaps(sequence, structure):
    remove = [i for i, nt in enumerate(sequence) if nt == '-']
    new_ss = [structure[i] for i in range(len(structure)) if i not in remove]
    return ''.join(new_ss)


def add_gaps(sequence, structure):
    structure = list(structure)
    new_ss = [structure.pop(0) if nt != '-' else '-' for  nt in sequence]
    return ''.join(new_ss)


def target_frequency(sequence, target):
    ss = remove_positioned_gaps(sequence, target)
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    return fc.pr_structure(ss)


def calculate_nts(sequence, region):
    return (region[1] - region[0]) - sequence[region[0]:region[1]].count('-')


def calculate_bps(sequence, region, symbol):
    return sequence[region[0]:region[1]].count(symbol)


def nt_in_structures(sequence, model_input):
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']
    # st1 = calculate_nts(sequence, stems[0][0][0]) + calculate_nts(sequence, stems[0][0][1]) + calculate_nts(sequence, stems[0][1][0]) + calculate_nts(sequence, stems[0][1][1])
    st2 = calculate_nts(sequence, stems[1][0]) + calculate_nts(sequence, stems[1][1])
    st3 = calculate_nts(sequence, stems[2][0]) + calculate_nts(sequence, stems[2][1])
    hl2 = calculate_nts(sequence, loops['hl2'])
    hl3 = calculate_nts(sequence, loops['hl3'])
    # upk1 = count_gaps(sequence, loops['upk1'])
    sI = calculate_nts(sequence, [stems[0][0][0][0],  stems[0][0][1][1]]) + calculate_nts(sequence, [stems[0][1][0][0],  stems[0][1][1][1]])
    sII = st2 + hl2
    sIII = st3 + hl3
    return  st2, st3, hl2, hl3, sI, sII, sIII # st1,


def bp_in_structures(sequence, model_input):
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']

    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    (ss, mfe) = fc.mfe()
    ss = add_gaps(sequence, ss)

    st1 = calculate_bps(ss, stems[0][0][0], '(') + calculate_bps(ss, stems[0][0][1], '(') + calculate_bps(ss, stems[0][1][0], ')') + calculate_bps(ss, stems[0][1][1], ')')
    st2 = calculate_bps(ss, stems[1][0], '(') + calculate_bps(ss, stems[1][1], ')')
    st3 = calculate_bps(ss, stems[2][0], '(') + calculate_bps(ss, stems[2][1], ')')
    hl2 = calculate_bps(ss, loops['hl2'], '.')
    hl3 = calculate_bps(ss, loops['hl3'], '.')

    sI = calculate_bps(ss, [stems[0][0][0][1],  stems[0][0][1][0]], '.') + calculate_bps(ss, [stems[0][1][0][1],  stems[0][1][1][0]], '.') + st1
    sII = st2 + hl2
    sIII = st3 + hl3
    return st1, st2, st3, hl2, hl3, sI, sII, sIII


def correct_bp_in_structure(ss, model_input, stems_pairs, predicted_pairs):
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']
    stem_counter = {0:0, 1:0, 2:0}
    for i, stem in enumerate(stems_pairs.values()):
        for pair in predicted_pairs:
            if pair in stem:
                stem_counter[i] += 1
    
    st1 = stem_counter[0] * 2
    st2 = stem_counter[1] * 2
    st3 = stem_counter[2] * 2

    hl2 = calculate_bps(ss, loops['hl2'], '.')
    hl3 = calculate_bps(ss, loops['hl3'], '.')

    sI = calculate_bps(ss, [stems[0][0][0][1],  stems[0][0][1][0]], '.') + calculate_bps(ss, [stems[0][1][0][1],  stems[0][1][1][0]], '.') + st1
    sII = st2 + hl2
    sIII = st3 + hl3
    return st1, st2, st3, hl2, hl3, sI, sII, sIII


def db_to_pairs(db):
    pairs = []
    db_dict = {')': '(',
               ']': '[',
               '}': '{',
               '>': '<',
               'a': 'A'}
    db_stacks = {'(':[],
               '[':[],
               '{':[],
               '<':[],
               'A': []
               }

    for i, symbol in enumerate(db):
        if symbol in db_stacks.keys():
            db_stacks[symbol].append(i)
        elif symbol in ')}]>a':
            opening_index = db_stacks[db_dict[symbol]].pop()
            closing_index = i
            pairs.append((opening_index,closing_index))
    if any(db_stacks.values()):
        print(db)
        print(db_stacks)
        raise ValueError(f"structure {db} not balanced ")
    return pairs


def get_pair(x, pairs):
    for pair in pairs:
        if x in pair:
            return pair


def get_stem_pairs(target_pairs, stems):
    stems_pairs = {}
    for stem_i, stem in enumerate(stems.values()):
        stem_x_pairs = []
        for stem_part in stem:
            if type(stem_part) == list:
                for j in stem_part:
                    for i in range(*j):
                        pair = get_pair(i, target_pairs)
                        stem_x_pairs.append(pair)
            else:
                for i in range(*stem_part):
                    pair = get_pair(i, target_pairs)
                    stem_x_pairs.append(pair)
        stem_x_pairs = list(set(stem_x_pairs))
        stems_pairs[stem_i] = stem_x_pairs
    return stems_pairs


def new_optimization_function(sequence, model_input, stems_pairs):
    target = model_input.structures[0]
    fc = RNA.fold_compound(sequence.replace('-',''))
    fc.pf()
    (ss, mfe) = fc.mfe()
    ss = add_gaps(sequence, ss)
    predicted_pairs = db_to_pairs(ss)

    st1, st2, st3, hl2, hl3, sI, sII, sIII = correct_bp_in_structure(ss, model_input, stems_pairs, predicted_pairs)

    cur_relations = {}
    cur_relations['st2_hl2'] = st2 / hl2
    cur_relations['st3_hl3'] = st3 / hl3

    cur_relations['sI_sII'] = sI / sII
    cur_relations['sI_sIII'] = sI / sIII
    cur_relations['sII_sIII'] = sII / sIII

    error_rel = 0
    error_rel = sum([error(relations[key][0], cur_relations[key]) for key in cur_relations.keys()])

    fc_pr = target_frequency(sequence, target)

    tot_length = len(sequence.replace('-', ''))
    error_len = error(target_len, tot_length)
    # relation_fraction_ = 1 / fc_pr
    new_val = fc_pr - (error_rel * fc_pr) - (error_len * length_fraction)

    # print('freq',fc_pr)
    # print('e rel',error_rel)
    # print('e len', error_len)
    # print('obj',new_val)
    return new_val, error_rel , error_len , fc_pr


def optimization_function(sequence, model_input):
    target = model_input.structures[0]
    st1, st2, st3, hl2, hl3, sI, sII, sIII = nt_in_structures(sequence, model_input)
    # st1, st2, st3, hl2, hl3, sI, sII, sIII = bp_in_structures(sequence, model_input)

    cur_relations = {}
    cur_relations['st2_hl2'] = st2 / hl2
    cur_relations['st3_hl3'] = st3 / hl3

    cur_relations['sI_sII'] = sI / sII
    cur_relations['sI_sIII'] = sI / sIII
    cur_relations['sII_sIII'] = sII / sIII

    error_rel = 0
    error_rel = sum([error(relations[key][0], cur_relations[key]) for key in cur_relations.keys()])

    fc_pr = target_frequency(sequence, target)

    tot_length = len(sequence.replace('-', ''))
    error_len = error(target_len, tot_length)
    new_val = max(0, fc_pr - (error_rel * relation_fraction) - (error_len * length_fraction))

    # print('freq',fc_pr)
    # print('e rel',error_rel)
    # print('e len', error_len)
    # print('obj',new_val)
    # return fc_pr
    return new_val, error_rel , error_len , fc_pr


def mc_optimize(model, model_input, objective, steps, temp, start=None):
    sampler = ir.Sampler(model)
    cur = sampler.sample() if start is None else start

    cur_seq = rna.values_to_seq(cur.values()[:len(model_input.structures[0])])
    curval, error_rel, error_len, fc_pr = objective(cur_seq)
    best, bestval = cur, curval
    ccs = model.connected_components()
    weights = [1/len(cc) for cc in ccs]

    values = []
    best_values = []
    errors_lens = []
    errors_rels = []
    freqs = []



    for i in tqdm((range(steps))): #tqdm
        cc = random.choices(ccs,weights)[0]
        new = sampler.resample(cc, cur)
        new_seq = rna.values_to_seq(new.values()[:len(model_input.structures[0])])
        newval, error_rel, error_len, fc_pr = objective(new_seq)

    
        if (newval >= curval or random.random() <= math.exp((newval-curval)/temp)):
            cur, curval = new, newval
            if curval > bestval:
                # print(curval)
                best, bestval = cur, curval
        if i%plotting_steps==0:
            values.append(newval)
            errors_lens.append(error_len)
            errors_rels.append(error_rel)
            freqs.append(fc_pr)
            best_values.append(bestval)

    # creating plot to have an overview of the behaviour of optimization function, length and relation error, frequency
    if False:
        x_ticks = [i * plotting_steps for i in range(len(values))]               
        plt.plot(x_ticks, values, label = "cur value", alpha=0.7) 
        plt.plot(x_ticks, best_values, label = "best value", alpha=0.7)
        plt.title('optimization function')
        plt.xlabel('steps')
        plt.legend()
        #plt.clf()
        plt.savefig('/scr/aldea/kgutenbrunner/working/xrRNA_design/TBFV_design/img/opt_function.png')

        plt.clf()
        plt.plot(x_ticks, values, alpha=0.7) 
        plt.plot(x_ticks, freqs, alpha=0.7)
        plt.plot(x_ticks, errors_rels, alpha=0.7)
        plt.plot(x_ticks, errors_lens, alpha=0.7)

        plt.title('total errors/score of the optimization function')
        plt.xlabel('steps')
        
        plt.legend(["obj function", "frequency", 'relation error', 'length error'])
        plt.savefig('/scr/aldea/kgutenbrunner/working/xrRNA_design/TBFV_design/img/errors.png')


        plt.clf()
        plt.plot(x_ticks, values, alpha=0.7) 
        plt.plot(x_ticks, freqs, alpha=0.7)
        plt.plot(x_ticks, np.array(errors_rels) * relation_fraction, alpha=0.7)
        plt.plot(x_ticks, np.array(errors_lens) * length_fraction, alpha=0.7)

        plt.title('relative errors/score of the optimization function')
        plt.xlabel('steps')
        
        plt.legend(["obj function", "frequency", f'relation error * {relation_fraction}', f'length error * {length_fraction}'])

        # plt.legend(["obj function", "frequency", f'relation error * {relation_fraction}', f'length error * {length_fraction}'])
        plt.savefig('/scr/aldea/kgutenbrunner/working/xrRNA_design/TBFV_design/img/errors_fraction.png')
    return (best, bestval), sampler


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
    sI_lens = []
    sII_lens = []
    sIII_lens = []

    for sample in samples:
        sample = rna.values_to_seq(sample.values()[:len(model_input.structures[0])])
        print(sample)
        print(model_input.structures[0])

        stems = model_input.structure_span['stem']
        loops = model_input.structure_span['loop']

        st1 = calculate_nts(sample, stems[0][0][0]) + calculate_nts(sample, stems[0][0][1]) + calculate_nts(sample, stems[0][1][0]) + calculate_nts(sample, stems[0][1][1])
        st2 = calculate_nts(sample, stems[1][0]) + calculate_nts(sample, stems[1][1])
        st3 = calculate_nts(sample, stems[2][0]) + calculate_nts(sample, stems[2][1])
        hl2 = calculate_nts(sample, loops['hl2'])
        hl3 = calculate_nts(sample, loops['hl3'])
        upk1 = calculate_nts(sample, loops['upk1'])
        
        sI = calculate_nts(sample, [stems[0][0][0][0],  stems[0][0][1][1]]) + calculate_nts(sample, [stems[0][1][0][0],  stems[0][1][1][1]])
        sII = st2 + hl2
        sIII = st3 + hl3

        st1_lens.append(st1)
        hl2_lens.append(hl2)
        st2_lens.append(st2)
        hl3_lens.append(hl3)
        st3_lens.append(st3)
        upk1_lens.append(upk1)
        sI_lens.append(sI)
        sII_lens.append(sII)
        sIII_lens.append(sIII)
        tot_lens.append(len(sample.replace('-', '')))

    print('-' * 150)
    print('ST1\t', sum(st1_lens) / len(st1_lens))
    print('HL2\t', sum(hl2_lens) / len(hl2_lens))
    print('ST2\t', sum(st2_lens) / len(st2_lens))
    print('HL3\t', sum(hl3_lens) / len(hl3_lens))
    print('ST3\t', sum(st3_lens) / len(st3_lens))
    print('s I\t', sum(sI_lens) / len(sI_lens))
    print('s II\t', sum(sII_lens) / len(sII_lens))
    print('s III\t', sum(sIII_lens) / len(sIII_lens))
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

