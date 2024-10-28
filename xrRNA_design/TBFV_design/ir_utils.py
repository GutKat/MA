import infrared as ir
import infrared.rna as rna
from collections import namedtuple
import math

ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_stem_regions var_loop_regions gaps target_length")

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
    'DiffInRange', 
    lambda pos_i, pos_j, A, B, var, name:  var([(name, pos_i), (name, pos_j)]),
    lambda y0, y1, A, B: 1 if A<= (y1-y0) and (y1-y0) <= B else 0,
    "ir_utils"
)
ir.def_constraint_class(
    'DiffInRangeAlpha',
    lambda pos_i, pos_j, pos_k, pos_l, A, B, var, name:  var([(name, pos_i), (name, pos_j), (name, pos_k), (name, pos_l)]),
    lambda y0, y1, y2, y3, A, B: 1 if A <= ((y1-y0) + (y3-y2)) <= B else 0, # A <=  ((y1-y0) + (y2-y3)) <= B
    "ir_utils"
)
########## END - Y for counting NT over possible gaps ###################


########## START - Y for calculating relations ###################
ir.def_constraint_class(
    # if two continous regions are calculated
    'CalculateRelationI',
    lambda pos_I, I_min, pos_II, II_min, relation, var, name:  var([(name, pos_I[0]), (name, pos_I[1]), (name, pos_II[0]), (name, pos_II[1])]),
    lambda y0, y1, x0, x1, I_min, II_min, relation: 1 if ((x1-x0) + II_min) and number_in_relation( (((y1-y0) + I_min) / ((x1-x0) + II_min)), relation) else 0,
    "ir_utils"
)
ir.def_constraint_class(
    # if one region is not continous (first one is seperated)
    'CalculateRelationII',
    lambda pos_I, I_min, pos_II, II_min, relation, var, name:  var([(name, pos_I[0]), (name, pos_I[1]), (name, pos_I[2]), (name, pos_I[3]), (name, pos_II[0]), (name, pos_II[1])]),
    lambda y0, y1, y2, y3, x0, x1, I_min, II_min, relation: 1 if ((x1-x0) + II_min) and number_in_relation( (((y1-y0) + (y3-y2) + I_min) / ((x1 - x0) + II_min)), relation) else 0,
    "ir_utils"
)
########## END - Y for calculating relations ###################



########## Start - Energy class ###################
ir.def_function_class(
    'BPEnergy',
    lambda i, j, is_terminal: [i, j],
    lambda x, y, is_terminal: bpenergy(x, y, is_terminal),
    "ir_utils"
)
########## End - Energy class ###################


def number_in_relation(x, relation):
    range_ = relations_range[relation]
    if range_[0] <= x <= range_[1]:
        return 1
    return 0


def bpenergy(x, y, is_terminal=False):
    if x == 4 or y == 4:
        return 0
    bpidx = bpindex_tab[x][y]
    return (params_bp[bpidx//2 + (3 if is_terminal else 0)]
            if bpidx >= 0 else -math.inf)


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

# min max
relations_range = {
        'stemG_hlG': [0.6, 1],
        'A_B': [2.0, 2.75],
        'A_G': [1.3, 1.4],
        'B_G': [0.5, 0.7]
    }



def create_model(model_input):
    ################## CREATE FEATURES ##################

    model_target = create_features(model_input)
    sampler = ir.Sampler(model_target)
    # min = 84, max = 90, mean = 86
    # need to add 2, have 2 nts in the beginning
    # 86 - 92
    sampler.set_target( 88, 4, 'totLength')
    sampler.set_target( -30, 10, 'energy')
    samples = [sampler.targeted_sample() for _ in range(1000)]
    feature_weights = sampler.model.features

    ################## MODEL SET-UP ##################

    n = len(model_input.structures[0])
    model = ir.Model()
    model.add_variables(n, 5, name='X')
    target_length = model_input.target_length
    iupac = str(model_input.iupac)

    bps = []
    for structure in model_input.structures:
        bps += rna.parse(structure)

    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()

    ################## CONSTRAINTS ##################

    ######## general constraints ########

    #add iupac of consensus sequence
    for i, x in enumerate(iupac):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))

    #add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)

    # make sure gaps in variable stems are only on the inner side --> to make sure we do not consider duplets
    for stem_a, stem_b in model_input.var_stem_regions:
        x, y = stem_a
        model.add_constraints(GapsRight(i) for i in range(x, y))

        x,y = stem_b
        model.add_constraints(GapsLeft(i) for i in range(x, y))

    #make sure gaps in variable loops are only on the right side --> to make sure we do not consider duplets
    for x, y in model_input.var_loop_regions:
        model.add_constraints(GapsRight(i) for i in range(x, y))

    # energy constraint
    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in bps], 'energy')
    model.set_feature_weight(feature_weights['energy'].weight, 'energy')

    ######## length constraints  ########

    ## target length with possible gaps
    n_y = len(possible_gaps)+1
    model.add_variables(n_y, n_y+1, name='Y')

    model.restrict_domains([('Y',0)], (0, 0))
    for i in range(1, n_y):
        model.restrict_domains([('Y',i)], (0, i-1))
    var = model.idx

    for i, g in enumerate(possible_gaps, start=1):
        model.add_constraints(IncrementYIfNotGap(g, i, var))


    gaps = model_input.gaps

    min_max_alpha = [2, 8]
    min_max_beta = [1, 6]
    min_max_gamma = [6,8]

    alpha_pos = [possible_gaps.index(gaps['alpha'][i]) if i%2==0 else possible_gaps.index(gaps['alpha'][i])+1 for i in range(len(gaps['alpha'])) ]
    beta_pos = [possible_gaps.index(gaps['beta'][i]) if i%2==0 else possible_gaps.index(gaps['beta'][i])+1 for i in range(len(gaps['beta'])) ]
    gamma_pos = [possible_gaps.index(gaps['gamma'][i]) if i%2==0 else possible_gaps.index(gaps['gamma'][i])+1 for i in range(len(gaps['gamma'])) ]
    stem_gamma_pos = [possible_gaps.index(gaps['gamma_stem'][i]) if i%2==0 else possible_gaps.index(gaps['gamma_stem'][i])+1 for i in range(len(gaps['gamma_stem'])) ]
    hl_gamma_pos = [possible_gaps.index(gaps['gamma_hl'][i]) if i%2==0 else possible_gaps.index(gaps['gamma_hl'][i])+1 for i in range(len(gaps['gamma_hl'])) ]

    # cur_min_alpha = 22
    # limit alpha to be between cur_min_alpha + 2 and cur_min_alpha + 8 --> 24 and 30
    model.add_constraints(DiffInRangeAlpha(*alpha_pos, *min_max_alpha, var, 'Y'))

    # cur_min_beta = 9
    # limit beta to be between cur_min_beta + 1 and cur_min_beta + 6 --> 10 and 15

    model.add_constraints(DiffInRange(*beta_pos, *min_max_beta, var, 'Y'))

    # min_gamma = 12
    # limit gamma to be between cur_min_gamma + 6 and cur_min_gamma + 8 --> 18 and 20
    model.add_constraints(DiffInRange(*gamma_pos, *min_max_gamma, var, 'Y'))

    # calculate relations between structure part
    model.add_constraints(CalculateRelationII(alpha_pos, 22, beta_pos, 9, 'A_B', var, 'Y')) # alpha / beta
    model.add_constraints(CalculateRelationII(alpha_pos, 22, gamma_pos, 12, 'A_G', var, 'Y')) # alpha / gamma
    model.add_constraints(CalculateRelationI(beta_pos, 9, gamma_pos, 12, 'B_G', var, 'Y')) # beta / gamma
    model.add_constraints(CalculateRelationII(stem_gamma_pos, 8, hl_gamma_pos, 4, 'stemG_hlG', var, 'Y')) # stem gamma / hl gamma

    if target_length:
        n_no_gaps = (len(possible_gaps) + target_length) - n
        model.add_constraints(CheckLastY(n_y, n_no_gaps, var))
    else:
        # Total Length
        model.add_functions([NotGap(i) for i in range(n)], 'totLength')
        model.set_feature_weight(feature_weights['totLength'].weight, 'totLength')

    return model

def create_features(model_input):

    ################## MODEL SET-UP ##################
    n = len(model_input.structures[0])
    model = ir.Model()
    model.add_variables(n, 5, name='X')
    target_length = model_input.target_length
    iupac = str(model_input.iupac)

    bps = []
    for structure in model_input.structures:
        bps += rna.parse(structure)

    possible_gaps = [i for i in range(len(iupac)) if iupac[i] == 'X']
    possible_gaps.sort()

    ################## CONSTRAINTS ##################
    ######## general constraints ########
    #add iupac of consensus sequence
    for i, x in enumerate(iupac):
        model.add_constraints(ir.ValueIn(i, iupacvalues(x)))

    #add basepair constraint of the structures
    for target in model_input.structures:
        ss = rna.parse(target)
        model.add_constraints(BPComp(i, j) for (i, j) in ss)

    # make sure gaps in variable stems are only on the inner side --> to make sure we do not consider duplets
    for stem_a, stem_b in model_input.var_stem_regions:
        x, y = stem_a
        model.add_constraints(GapsRight(i) for i in range(x, y))

        x,y = stem_b
        model.add_constraints(GapsLeft(i) for i in range(x, y))

    #make sure gaps in variable loops are only on the right side --> to make sure we do not consider duplets
    for x, y in model_input.var_loop_regions:
        model.add_constraints(GapsRight(i) for i in range(x, y))

    # energy constraint
    model.add_functions([BPEnergy(i, j, (i-1, j+1) not in bps) for (i,j) in bps], 'energy')

    ######## length constraints  ########

    ## target length with possible gaps
    n_y = len(possible_gaps)+1
    model.add_variables(n_y, n_y+1, name='Y')

    model.restrict_domains([('Y',0)], (0, 0))
    for i in range(1, n_y):
        model.restrict_domains([('Y',i)], (0, i-1))
    var = model.idx

    for i, g in enumerate(possible_gaps, start=1):
        model.add_constraints(IncrementYIfNotGap(g, i, var))

    gaps = model_input.gaps

    min_max_alpha = [2, 8]
    min_max_beta = [1, 6]
    min_max_gamma = [6,8]

    alpha_pos = [possible_gaps.index(gaps['alpha'][i]) if i%2==0 else possible_gaps.index(gaps['alpha'][i])+1 for i in range(len(gaps['alpha'])) ]
    beta_pos = [possible_gaps.index(gaps['beta'][i]) if i%2==0 else possible_gaps.index(gaps['beta'][i])+1 for i in range(len(gaps['beta'])) ]
    gamma_pos = [possible_gaps.index(gaps['gamma'][i]) if i%2==0 else possible_gaps.index(gaps['gamma'][i])+1 for i in range(len(gaps['gamma'])) ]
    stem_gamma_pos = [possible_gaps.index(gaps['gamma_stem'][i]) if i%2==0 else possible_gaps.index(gaps['gamma_stem'][i])+1 for i in range(len(gaps['gamma_stem'])) ]
    hl_gamma_pos = [possible_gaps.index(gaps['gamma_hl'][i]) if i%2==0 else possible_gaps.index(gaps['gamma_hl'][i])+1 for i in range(len(gaps['gamma_hl'])) ]

    # cur_min_alpha = 22
    # limit alpha to be between cur_min_alpha + 2 and cur_min_alpha + 8 --> 24 and 30
    model.add_constraints(DiffInRangeAlpha(*alpha_pos, *min_max_alpha, var, 'Y'))

    # cur_min_beta = 9
    # limit beta to be between cur_min_beta + 1 and cur_min_beta + 6 --> 10 and 15
    model.add_constraints(DiffInRange(*beta_pos, *min_max_beta, var, 'Y'))

    # min_gamma = 12
    # limit gamma to be between cur_min_gamma + 6 and cur_min_gamma + 8 --> 18 and 20
    model.add_constraints(DiffInRange(*gamma_pos, *min_max_gamma, var, 'Y'))

    # calculate relations between structure part
    model.add_constraints(CalculateRelationII(alpha_pos, 22, beta_pos, 9, 'A_B', var, 'Y')) # alpha / beta
    model.add_constraints(CalculateRelationII(alpha_pos, 22, gamma_pos, 12, 'A_G', var, 'Y')) # alpha / gamma
    model.add_constraints(CalculateRelationI(beta_pos, 9, gamma_pos, 12, 'B_G', var, 'Y')) # beta / gamma
    model.add_constraints(CalculateRelationII(stem_gamma_pos, 8, hl_gamma_pos, 4, 'stemG_hlG', var, 'Y')) # stem gamma / hl gamma



    if target_length:
        n_no_gaps = (len(possible_gaps) + target_length) - n
        model.add_constraints(CheckLastY(n_y, n_no_gaps, var))
    else:
        # Total Length
        model.add_functions([NotGap(i) for i in range(n)], 'totLength')

    return model

