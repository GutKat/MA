import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import itertools as it
from collections import Counter


ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_stem_regions var_loop_regions structure_span target_length")


relations_range = {
        'st2_hl2': [2, 4],
        'st3_hl3': [3/4, 2],
        'sI_sII': [1.5, 2.6],
        'sII_sIII': [0.5, 0.9],
        'sI_sIII': [1.20, 1.4]
    }

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

ir.def_constraint_class( 
    'IncrementIfBP',
    lambda i, j, var: var([('Y',i),('Y',j)]),
    lambda x, y: print(x,y), #z+1 if (x,y) in _bpcomp_tab else
    "ir_utils"
) 

ir.def_function_class( 
    'NotGap',
    lambda i: [i],
    lambda x: 1 if x != 4 else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'StartYIfNotGap',
    lambda i, start, var: [*var([('X',i)]), *var([('Y',start)])],
    lambda x, y, start: 1 if (start==y and x==4) or (start+1==y and x!=4) else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'CheckLastY',
    lambda i, A, var: var([('Y',i-1)]),
    lambda y, A: 1 if y == A else 0,# p1 if x <= cutOff else 0 print(x, c)
    "ir_utils"
)


ir.def_constraint_class(
    'IncrementYIfNotGap',
    lambda x, pos, var: [*var([('X',x)]), *var([('Y',pos-1), ('Y',pos)])],
    lambda nt, y0, y1: 1 if (nt != 4 and y0 + 1 == y1) or (nt == 4 and y0 == y1) else 0,
    "ir_utils"
)


def number_in_range(x, range_):
    if range_[0] <= round(x, 3) <= range_[1]:
        return 1
    return 0


def create_model(model_input):


    n = len(model_input.structures[0])

    model = ir.Model()
    model.add_variables(n,5)
    target_length = model_input.target_length
    iupac = str(model_input.iupac)

    #add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)
        #model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in ss)], 'energy')

    #add iupac of consensus sequence
    for i, x in enumerate(iupac):
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


    # new variable for the target length constraint

    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()
    n_y = len(possible_gaps)
    model.add_variables(n_y, n_y+1, name='Y')
    for i in range(n_y):
        model.restrict_domains([('Y',i)], (0, i))
    var = model.idx

    model.add_constraints(StartYIfNotGap(possible_gaps[0], 0, var))

    for i, g in enumerate(possible_gaps[1:], start=1):
        model.add_constraints(IncrementYIfNotGap(g, i, var))

    n_no_gaps = (len(possible_gaps) + target_length) - n
    model.add_constraints(CheckLastY(n_y, n_no_gaps, var))
 
    # ######## constraints for relation controll ########
    # stems = model_input.structure_span['stem']
    # loops = model_input.structure_span['loop']
    # print(stems)
    # print(loops)
    # # check helix beta and hairpin beta relation
    # hlII = loops['hl2']
    # stII = stems[1]
    # model.add_constraints(Testing(hlII, stII, var))

    # ######## additional constraints from simulation analysis ########

    # loops = model_input.structure_span['loop']
    # hl2 = range(*loops['hl2'])
    # upk1 = range(*loops['upk1'])

    return model
