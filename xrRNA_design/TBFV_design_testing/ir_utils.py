import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import itertools as it
from collections import Counter


ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_stem_regions var_loop_regions structure_span target_length")


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

# pairs     = [(A, U), (C, G), (G, C), (G, U), (U, A), (U, G)]
_bpcomp_tab = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
_bpcomp_tab_gap = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2), (4,4)]


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

ir.def_constraint_class(
    'StartZIfBP',
    lambda x_i, x_j, start, var: [*var([('X',x_i), ('X',x_j)]), *var([('Z',start)])], #[*var([('X',x_i)]), *var([('X',x_j)]), *var([('Z',z_pos-1), ('Z',z_pos)])],
    lambda ntA, ntB, z0, start: 1 if ((ntA, ntB) != (0,3) and z0 == start) or ((ntA, ntB) == (0,3) and z0 == start + 1) else 0, #print(ntA, ntB, z0, z1) if ((ntA, ntB) not in _bpcomp_tab and z0 == z1) or ((ntA, ntB) in _bpcomp_tab and z0 + 1 == z1) else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'IncrementZIfBP',
    lambda x_i, x_j, z_pos, var: [*var([('X',x_i), ('X',x_j)]), *var([('Z',z_pos-1), ('Z',z_pos)])], #[*var([('X',x_i)]), *var([('X',x_j)]), *var([('Z',z_pos-1), ('Z',z_pos)])],
    lambda ntA, ntB, z0, z1: 1 if ((ntA, ntB) not in _bpcomp_tab and z0 == z1) or ((ntA, ntB) in _bpcomp_tab and z0 + 1 == z1) else 0, #print(ntA, ntB, z0, z1) if ((ntA, ntB) not in _bpcomp_tab and z0 == z1) or ((ntA, ntB) in _bpcomp_tab and z0 + 1 == z1) else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'CheckLastZ',
    lambda i, A, var: var([('Z',i-1)]),
    lambda z, A: 1 if z == A else 0,# p1 if x <= cutOff else 0 print(x, c)
    "ir_utils"
)

def create_model(model_input):


    n = len(model_input.structures[0])
    model = ir.Model()
    model.add_variables(n, 5, name='X')
    target_length = model_input.target_length
    iupac = str(model_input.iupac)

    #add iupac of consensus sequence
    for i, x in enumerate(iupac):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))


    ######## additional constraints from simulation analysis ########

    # target length as network function
    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()
    model.add_functions([NotGap(i) for i in possible_gaps], 'gaps')




    # if 'X' in iupac:
    #     possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    #     possible_gaps.sort()
    #     n_y = len(possible_gaps)
    #     model.add_variables(n_y, n_y+1, name='Y')
    #     var = model.idx
    #     print('possible gaps:', n_y)
    #
    #
    #     model.add_constraints(StartYIfNotGap(possible_gaps[0], 0, var))
    #
    #     for i, g in enumerate(possible_gaps[1:], start=1):
    #         model.add_constraints(IncrementYIfNotGap(g, i, var))
    #
    #     n_no_gaps = (len(possible_gaps) + target_length) - n
    #     print('target gaps:', n - target_length)
    #     print('target no gaps:', n_no_gaps)
    #     model.add_constraints(CheckLastY(n_y, n_no_gaps, var))


    # adding bp constraints

    # loops = model_input.structure_span['loop']
    # hl2 = range(*loops['hl2'])
    # upk1 = range(*loops['upk1'])
    # print(hl2, upk1)


    #
    # combs = []
    # for pos_i in hl2:
    #     for pos_j in upk1:
    #         combs.append((pos_i, pos_j))
    # n_z = len(combs)
    # model.add_variables(n_z, n_z+1, name='Z')
    # var = model.idx
    #
    # model.add_constraints(StartZIfBP(*combs[0], 0, var))
    # for i, (pos_i, pos_j) in enumerate(combs[1:], start=1):
    #     print(pos_i, pos_j, i)
    #     model.add_constraints(IncrementZIfBP(pos_i, pos_j, i, var))
    #
    # AU_bp = 3
    # model.add_constraints(CheckLastZ(n_z, AU_bp, var))
    # print('model')


    return model
