import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import math

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

params_bp = [
    -0.52309,
    -2.10208,
    -0.88474,
    1.2663,
    -0.0907,
    0.78566
]

params_stacking = [
    -0.18826,
    -0.2651,
    -1.13291,
    -1.09787,
    -0.38606,
    -0.62086,
    -1.11752,
    -1.10548,
    -2.2374,
    -1.89434,
    -1.22942,
    -1.44085,
    -0.55066,
    -0.49625,
    -1.26209,
    -1.58478,
    -0.72185,
    -0.68876
]

bpindex_tab = [[-1, -1, -1, 0],
               [-1, -1, 2, -1],
               [-1, 3, -1, 4],
               [1, -1, 5, -1]]

# pairs     = [(A, U), (C, G), (G, C), (G, U), (U, A), (U, G)]
_bpcomp_tab = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
_bpcomp_tab_gap = [(0, 3), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2), (4,4)]


def iupacvalues(symbol):
    return [ rna.nucleotide_to_value(x) for x in extended_iupac_nucleotides[symbol] ]


def bpenergy(x, y, is_terminal=False):
    if x == 4 or y == 4:
        return 0
    bpidx = bpindex_tab[x][y]
    return (params_bp[bpidx//2 + (3 if is_terminal else 0)]
            if bpidx >= 0 else -math.inf)


########## Start - BP/NotBP/Gap/NotGap ###################

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
    'Gap',
    lambda i: [i],
    lambda x: 1 if x == 4 else 0,
    "ir_utils"
)

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

########## End - BP/NotBP/Gap/NotGap ###################


########## START - Y for counting NT over possible gaps ###################

ir.def_constraint_class(
    'Is',
    lambda i, A, var, name: var([(name,i-1)]),
    lambda y, A: 1 if y == A else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'InRange',
    lambda i, A, B, var, name: var([(name,i-1)]),
    lambda y, A, B: 1 if A <= y <= B else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'IncrementIfNotGap',
    lambda x, pos, var, name: [*var([('X',x)]), *var([(name,pos-1), (name,pos)])],
    lambda nt, y0, y1: 1 if (nt != 4 and y0 + 1 == y1) or (nt == 4 and y0 == y1) else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'SmallerThan',
    lambda i, A, var, name: var([(name,i-1)]),
    lambda z, A: 1 if z <= A else 0,
    "ir_utils"
)

ir.def_constraint_class(
    'DiffInRange', 
    lambda pos_i, pos_j, A, B, var, name:  var([(name, pos_i), (name, pos_j)]),
    lambda y0, y1, A, B: 1 if A<= (y1-y0) <= B else 0,
    "ir_utils"
)

########## END - Y for counting NT over possible gaps ###################



########## Start - Energy class ###################

ir.def_function_class(
    'BPEnergy',
    lambda i, j, is_terminal: [i, j],
    lambda x, y, is_terminal: bpenergy(x, y, is_terminal),
    "ir_utils"
)

########## End - Energy class ###################



def create_model(model_input, feature_weights):
    '''
    Args:
        model_input (dict): 
        feature_weights (dict): 

    Returns:
        model (infrared.model):
    '''
    ################## MODEL SET-UP ##################

    n = len(model_input.structures[0])
    model = ir.Model()
    model.add_variables(n, 5, name='X')
    target_length = model_input.target_length
    iupac = str(model_input.iupac)

    # get all base pairs from the input structures
    bps = []
    for structure in model_input.structures:
        bps += rna.parse(structure)

    # get all possible gap position for variable Y
    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()

    ################## CONSTRAINTS ##################

    ######## general constraints ########

    # add iupac of consensus sequence
    for i, x in enumerate(iupac):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))

    # add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)

    # make sure gaps in variable stems are only on the inner side --> to make sure we do not consider duplets
    for stem_a, stem_b in model_input.var_stem_regions:
        x, y = stem_a
        model.add_constraints(GapsRight(i) for i in range(x, y))

        x,y = stem_b
        model.add_constraints(GapsLeft(i) for i in range(x, y))


    # make sure gaps in variable loops are only on the right side --> to make sure we do not consider duplets
    for x, y in model_input.var_loop_regions:
        model.add_constraints(GapsRight(i) for i in range(x, y))

    # energy constraint
    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in bps], 'energy')
    model.set_feature_weight(feature_weights['energy'].weight, 'energy')

    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in rna.parse(model_input.structures[2])], 'energy_pk2')
    model.set_feature_weight(-1.5, 'energy_pk2')
    # model.set_feature_weight(feature_weights['energy_pk2'].weight, 'energy_pk2')


    ######## length constraints and constraints for hair pin loop beta ########

    ## target length with possible gaps
    min_beta = 2
    max_beta = 10
    min_gamma = 8
    max_gamma = 11
    min_hl_gamma_nt = 1
    max_hl_gamma_nt = 11
    hl_gamma = model_input.structure_span['loop']['hl3']

    if target_length:
        ss_beta = model_input.structure_span['gaps']['beta']
        ss_gamma = model_input.structure_span['gaps']['gamma']
        
        # set up variable for counting over possible gap position
        n_y = len(possible_gaps)+1
        model.add_variables(n_y, n_y+1, name='Y')

        # restrict Y for efficiency
        model.restrict_domains([('Y',0)], (0, 0))
        for i in range(1, n_y):
            model.restrict_domains([('Y',i)], (0, i-1))
        var = model.idx

        # go over possible gap position and increase Y if position is not a gap
        for i, pos in enumerate(possible_gaps, start=1):
            model.add_constraints(IncrementIfNotGap(pos, i, var, 'Y'))

        if type(target_length) == int:      
        # if target length should be exact number
            # check if Y is same as n_no_gaps
            n_no_gaps = (len(possible_gaps) + target_length) - n
            model.add_constraints(Is(n_y, n_no_gaps, var, 'Y'))
        
        elif type(target_length) == list:   
        # if target length should be within a range
            # check if Y is within the target length range (target_A, target_B)
            target_A = (len(possible_gaps) + target_length[0]) - n
            target_B = (len(possible_gaps) + target_length[1]) - n
            model.add_constraints(InRange(n_y, target_A, target_B, var, 'Y'))

        # limit beta to be between min_beta + 5 and max_beta + 5 --> 7 and 15
        pos_A = possible_gaps.index(ss_beta[0])         # plus one because we have an artificial zero in variable Y
        pos_B = possible_gaps.index(ss_beta[1]) + 1     # plus one because we have an artificial zero in variable Y
        model.add_constraints(DiffInRange(pos_A, pos_B, min_beta, max_beta, var, 'Y'))

        # limit gamma to be between min_gamma + 12 and max_gamma + 12 --> 20 and 23
        pos_A = possible_gaps.index(ss_gamma[0])         # plus one because we have an artificial zero in variable Y
        pos_B = possible_gaps.index(ss_gamma[1]) + 1     # plus one because we have an artificial zero in variable Y
        model.add_constraints(DiffInRange(pos_A, pos_B, min_gamma, max_gamma, var, 'Y'))

        # limit hairpin loop gamma to be between min_hl_gamma_nt + 4 and max_hl_gamma_nt + 4 --> 5 and 15
        pos_A = possible_gaps.index(hl_gamma[0])         # plus one because we have an artificial zero in variable Y
        pos_B = possible_gaps.index(hl_gamma[1]) + 1     # plus one because we have an artificial zero in variable Y
        model.add_constraints(DiffInRange(pos_A, pos_B, min_hl_gamma_nt, max_hl_gamma_nt, var, 'Y'))

    # if no target length was set, Y is not needed for counting gaps
    else:
        ss_beta = model_input.structure_span['ss']['beta']
        ss_gamma = model_input.structure_span['ss']['gamma']

        # get all possible gap position within beta and set up variable Z_B
        gaps_in_beta = [i for i in range(ss_beta[0],ss_beta[1]) if iupac[i] == 'X']
        n_beta = len(gaps_in_beta) + 1
        model.add_variables(n_beta, n_beta + 1, name='Z_B')
        model.restrict_domains([('Z_B', 0)], (0, 0))
        var = model.idx
        for i, pos in enumerate(gaps_in_beta, start=1):
            model.add_constraints(IncrementIfNotGap(pos, i, var, 'Z_B'))           
        # limit beta to be between min_beta + 5 and max_beta + 5 --> 7 and 15
        model.add_constraints(InRange(n_beta, min_beta, max_beta, var, 'Z_B'))

        # get all possible gap position within gamma and set up variable Z_G
        gaps_in_gamma = [i for i in range(ss_gamma[0],ss_gamma[1]) if iupac[i] == 'X']
        n_gamma = len(gaps_in_gamma) + 1
        model.add_variables(n_gamma, n_gamma + 1, name='Z_G')
        model.restrict_domains([('Z_G', 0)], (0, 0))
        var = model.idx
        for i, pos in enumerate(gaps_in_gamma, start=1):
            model.add_constraints(IncrementIfNotGap(pos, i, var, 'Z_G'))
        # limit hairpin loop gamma to be between min_gamma + 12 and max_gamma + 12 --> 20 and 23
        model.add_constraints(InRange(n_gamma, min_gamma, max_gamma, var, 'Z_G'))

        # limit hairpin loop gamma to be between min_hl_gamma_nt + 4 and max_hl_gamma_nt + 4 --> 5 and 15
        pos_A = gaps_in_gamma.index(hl_gamma[0])         # plus one because we have an artificial zero in variable Y
        pos_B = gaps_in_gamma.index(hl_gamma[1]) + 1     # plus one because we have an artificial zero in variable Y
        model.add_constraints(DiffInRange(pos_A, pos_B, min_hl_gamma_nt, max_hl_gamma_nt, var, 'Z_G'))

        # set weight for total Length
        model.add_functions([NotGap(i) for i in range(n)], 'totLength')
        model.set_feature_weight(feature_weights['totLength'].weight, 'totLength')

    return model


def create_model_target(model_input):
    '''
    same as function create model, but instead it sets no feature weights
    used to find weight for feature by target sampling
    for detailed information take a look at function "create_model"
    '''
    
    ################## MODEL SET-UP ##################
    n = len(model_input.structures[0])
    model = ir.Model()
    model.add_variables(n, 5, name='X')
    iupac = str(model_input.iupac)

    bps = []
    for structure in model_input.structures:
        bps += rna.parse(structure)

    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()

    ################## CONSTRAINTS ##################
    ######## general constraints ########
    # sequence constraints
    for i, x in enumerate(iupac):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))

    # base pairs constraints
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)

    # redundancy constraints
    for stem_a, stem_b in model_input.var_stem_regions:
        x, y = stem_a
        model.add_constraints(GapsRight(i) for i in range(x, y))

        x,y = stem_b
        model.add_constraints(GapsLeft(i) for i in range(x, y))

    for x, y in model_input.var_loop_regions:
        model.add_constraints(GapsRight(i) for i in range(x, y))

    # energy
    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in bps], 'energy')
    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in rna.parse(model_input.structures[2])], 'energy_pk2')

    ######## length constraints and constraints for hair pin loop beta ########
    min_beta = 2
    max_beta = 10
    min_gamma = 8
    max_gamma = 11
    min_hl_gamma_nt = 1
    max_hl_gamma_nt = 11
    hl_gamma = model_input.structure_span['loop']['hl3']

    ss_beta = model_input.structure_span['ss']['beta']
    ss_gamma = model_input.structure_span['ss']['gamma']

    # beta length
    gaps_in_beta = [i for i in range(ss_beta[0],ss_beta[1]) if iupac[i] == 'X']
    n_beta = len(gaps_in_beta) + 1
    model.add_variables(n_beta, n_beta + 1, name='Z_B')
    model.restrict_domains([('Z_B', 0)], (0, 0))
    var = model.idx
    for i, pos in enumerate(gaps_in_beta, start=1):
        model.add_constraints(IncrementIfNotGap(pos, i, var, 'Z_B'))           
    model.add_constraints(InRange(n_beta, min_beta, max_beta, var, 'Z_B'))

    # gamma length
    gaps_in_gamma = [i for i in range(ss_gamma[0],ss_gamma[1]) if iupac[i] == 'X']
    n_gamma = len(gaps_in_gamma) + 1
    model.add_variables(n_gamma, n_gamma + 1, name='Z_G')
    model.restrict_domains([('Z_G', 0)], (0, 0))
    var = model.idx
    for i, pos in enumerate(gaps_in_gamma, start=1):
        model.add_constraints(IncrementIfNotGap(pos, i, var, 'Z_G'))           
    model.add_constraints(InRange(n_gamma, min_gamma, max_gamma, var, 'Z_G'))

    # gamma hairpin length
    pos_A = gaps_in_gamma.index(hl_gamma[0]) + 1
    pos_B = gaps_in_gamma.index(hl_gamma[1]) + 1
    model.add_constraints(DiffInRange(pos_A, pos_B, min_hl_gamma_nt, max_hl_gamma_nt, var, 'Z_G'))

    # total Length
    model.add_functions([NotGap(i) for i in range(n)], 'totLength')
    return model

