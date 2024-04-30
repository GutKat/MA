import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import itertools as it


ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_stem_regions var_loop_regions structure_span")


extended_iupac_nucleotides = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'U',
    'U': 'U',
    'R': 'AG',
    'Y': 'CU',
    'S': 'CG',
    'W': 'AU',
    'K': 'GU',
    'M': 'AC',
    'B': 'CGU',
    'D': 'AGU',
    'H': 'ACU',
    'V': 'ACG',
    'N': 'ACGU',
    'X': 'ACGU-',
    '-': '-',
    '.': '-'
}


def iupacvalues(symbol):
    return [ rna.nucleotide_to_value(x) for x in extended_iupac_nucleotides[symbol] ]

_bpcomp_tab = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
_bpcomp_tab_gap = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2), (4,4)]

ir.def_constraint_class( 
    'GapsRight',
    lambda i: [i,i+1],
    lambda left,right: 0 if left == 4 and right != 4 else 1,
    "ir_utils"
    # NNN-- True
    # --NNN False
)


ir.def_constraint_class(
    'GapsLeft',
    lambda i: [i,i+1],
    lambda left,right: 0 if left != 4 and right == 4 else 1,
    "ir_utils"
    # --NNN True
    # NNN-- False
)


ir.def_function_class( 
    'NotGap',
    lambda i: [i],
    lambda x: 1 if x != 4 else 0,
    "ir_utils"

) 

ir.def_constraint_class( 
    'BPComp',
    lambda i,j: [i,j],
    lambda x,y: 1 if (x,y) in _bpcomp_tab_gap else 0,
    "ir_utils"
) 

ir.def_constraint_class( 
    'BPNotComp',
    lambda i,j: [i,j],
    lambda x,y: 0 if (x,y) in _bpcomp_tab_gap else 1,
    "ir_utils"
) 


def create_model(model_input):


    n = len(model_input.structures[0])

    model = ir.Model()
    model.add_variables(n,5)

    #add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)
        #model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in ss)], 'energy')

    #add iupac of consensus sequence
    for i, x in enumerate(str(model_input.iupac)):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))
    
    # make sure gaps in variable stems are only on the inner side --> to make sure we do not consider duplets
    for stem_a, stem_b in model_input.var_stem_regions:
        x, y = stem_a
        model.add_constraints(GapsRight(i) for i in range(x, y))

        x,y = stem_b
        model.add_constraints(GapsLeft(i) for i in range(x, y))

    #make sure gaps in variable loops are only on the right side --> to make sure we do not consider duplets
    for x, y in model_input.var_loop_regions:
        model.add_constraints(GapsRight(i) for i in range(x, y))


    ######## Manually added length-weights of the structures ########
    # stems = model_input.structure_span['stem']
    # loops = model_input.structure_span['loop']

    # # Stem I
    # model.add_functions([NotGap(i) for i in it.chain(range(*stems[0][0]), range(*stems[0][1]))], 'st1')
    # # Hairpin-loop II
    # model.add_functions([NotGap(i) for i in range(*loops['hl2'])], 'hl2')
    # # Stem II
    # model.add_functions([NotGap(i) for i in it.chain(range(*stems[1][0]), range(*stems[1][1]))], 'st2')
    # # Hairpin-loop III
    # model.add_functions([NotGap(i) for i in range(*loops['hl3'])], 'hl3')
    # # Stem III
    # model.add_functions([NotGap(i) for i in it.chain(range(*stems[2][0]), range(*stems[2][1]))], 'st3')
    # # nt until PK1
    # model.add_functions([NotGap(i) for i in range(*loops['upk1'])], 'upk1')


    # model.set_feature_weight(-0.5, 'st1')
    # model.set_feature_weight(4, 'hl2')
    # model.set_feature_weight(-0.5, 'st2')
    # model.set_feature_weight(-2, 'hl3')
    # model.set_feature_weight(-0.5, 'st3')
    # model.set_feature_weight(2, 'upk1')

    #Total Length of sequence
    model.add_functions([NotGap(i) for i in range(n)], 'totLength')

    return model

'''Stuff for testing the constraints'''

# testing
def gapslefttest(seq, region):
    embedding = {'N':0, '-': 4, 'X': 4}
    seq_embedded = [embedding[x] for x in seq]
    results = []
    x, y = region
    for i in range(x, y):
        left, right = seq_embedded[i], seq_embedded[i+1]
        results.append(0 if left != 4 and right == 4 else 1)
    return 0 if 0 in results else 1


def gapsrighttest(seq, region):
    embedding = {'N':0, '-': 4, 'X': 4}
    seq_embedded = [embedding[x] for x in seq]
    results = []
    x,y = region
    for i in range(x, y):
        left, right = seq_embedded[i], seq_embedded[i+1]
        results.append(0 if left == 4 and right != 4 else 1)
    return 0 if 0 in results else 1


def gapstemtest(seq, region_a, region_b):
    results = [gapsrighttest(seq, region_a), gapslefttest(seq, region_b)]
    return 0 if 0 in results else 1


def constraint_testing(test='test1'):
    if test == 'test1':
        test_seqs_loop = ['NN---', 'N---N', '--N--', '---NN']
        true_left_loop = [0, 0, 0, 1]
        true_right_loop = [1, 0, 0, 0]
        for i, seq in enumerate(test_seqs_loop):
            assert gapslefttest(seq, 0, 4) == true_left_loop[i]
            assert gapsrighttest(seq, 0, 4) == true_right_loop[i]


        test_seqs_stem = ['NN----NN', '--NNNN--']
        true_stem = [1, 0]
        for i, seq in enumerate(test_seqs_stem):
            assert gapstemtest(seq, (0, 3), (4,7)) == true_stem[i]

    else:
        seq, var_loops, var_stems = test
        results = []
        for loop in var_loops:
            results.append(gapsrighttest(seq, loop))
        for stem in var_stems:
            results.append(gapstemtest(seq, stem[0], stem[1]))
        return False if 0 in results else True



def main():
    constraint_testing()

if __name__ == "__main__":
    main()