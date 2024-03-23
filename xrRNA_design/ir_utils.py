import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import itertools as it


ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_stem_regions var_loop_regions")


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


ir.def_constraint_class( 
    'EmptyIfLeftEmpty',
    lambda i: [i,i+1],
    lambda left,right: 0 if left == 4 and right != 4 else 1,
    "ir_utils"
)

ir.def_constraint_class( 
    'BPorGap',
    lambda i,j: [i,j],
    lambda x,y: 1 if rna.BPComp(x,y) else 1 if (x == 4 and y == 4) else 0,
    "ir_utils"
)
 

ir.def_function_class( 
    'NotGap',
    lambda i: [i],
    lambda x: 1 if x != 4 else 0,
    "ir_utils"
) 


def create_model(model_input):


    n = len(model_input.structures[0])

    model = ir.Model()


    #First n variables descibe the base positions. Have 4 different bases and a gap (domain of five)
    model.add_variables(n,5)

    #Variables for the variable region.
    model.add_variables(n,n)

    #add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(rna.BPComp(i, j) for (i, j) in ss)
        #model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in ss)], 'energy')

    #add iupac of consensus sequence
    for i, x in enumerate(str(model_input.iupac)):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))
    
    #add the variable region of the stems
    for stem_a, stem_b in model_input.var_stem_regions:
        for i,j in zip(range(stem_a[0],stem_a[1]+1), range(stem_b[1],stem_b[0]-1, -1)):
            model.add_constraints(BPorGap(i, j))

    #add possible unpaired nt of the loops
    for x, y in model_input.var_loop_regions:
        model.add_constraints(EmptyIfLeftEmpty(i) for i in range(x, y))



    ######## Manually added length-weights of the structures ########
    stems = {0: [(0,16), (63, 79)], 1: [(22,27), (31, 36)], 2: [(37, 42), (57, 62)]}
    loops = {'hl2': (27,31), 'hl3': (42, 57), 'upk1': (79, 94)}

    # Stem I
    model.add_functions([NotGap(i) for i in it.chain(range(*stems[0][0]), range(*stems[0][1]))], 'st1')

    # Hairpin-loop II
    model.add_functions([NotGap(i) for i in range(*loops['hl2'])], 'hl2')
    # Stem II
    model.add_functions([NotGap(i) for i in it.chain(range(*stems[1][0]), range(*stems[1][1]))], 'st2')

    # Hairpin-loop III
    model.add_functions([NotGap(i) for i in range(*loops['hl3'])], 'hl3')
    # Stem III
    model.add_functions([NotGap(i) for i in it.chain(range(*stems[2][0]), range(*stems[2][1]))], 'st3')

    # nt until PK1
    model.add_functions([NotGap(i) for i in range(*loops['upk1'])], 'upk1')

    model.set_feature_weight(-3, 'st1')
    model.set_feature_weight(4, 'hl2')
    model.set_feature_weight(-1, 'st2')
    model.set_feature_weight(-3, 'hl3')
    model.set_feature_weight(-1, 'st3')
    model.set_feature_weight(3, 'upk1')

    #Total Length of sequence
    model.add_functions([NotGap(i) for i in range(n)], 'totLength')

    return model



