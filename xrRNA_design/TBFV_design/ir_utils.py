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
########## END - Y for counting NT over possible gaps ###################


########## START - Y for calculating relations ###################
ir.def_constraint_class(
    'CalculateRelationStemHairpin',
    lambda stem, stem_c, hl, hl_c, relation, var: [*var([('Y',stem[0]), ('Y',stem[1]), ('Y',stem[2]), ('Y',stem[3])]), *var([('Y',hl[0]), ('Y',hl[1])])],
    lambda x0, x1, x2, x3, y0, y1, relation, stem_c, hl_c: 1 if ((y0 - y1) + hl_c) and number_in_relation(((x0-x1)+(x2-x3)+stem_c*2) / ((y0 - y1) + hl_c), relation) else 0,
    "ir_utils"
)
ir.def_constraint_class(
    'CalculateRelationStructure',
    lambda stI, stI_c, stJ, stII_c, relation, var: [*var([('Y',stI[0]), ('Y',stI[1]), ('Y',stJ[0]), ('Y',stJ[1])])],
    lambda x0, x1, y0, y1, relation, stI_c, stII_c: 1 if ((y0 - y1) + stII_c) and number_in_relation(((x0-x1)+stI_c) / ((y0 - y1) + stII_c), relation) else 0,
    "ir_utils"
)
ir.def_constraint_class(
    'CalculateRelationStructureI',
    lambda stI, stI_c, stII, stII_c, relation, var: [*var([('Y',stI[0]), ('Y',stI[1]), ('Y',stI[2]), ('Y',stII[0]), ('Y',stII[1])])],
    lambda x0, x1, x2, y0, y1, relation, stI_c, stII_c: 1 if ((y0 - y1) + stII_c) and number_in_relation(((x0-x1)+(x2)+stI_c*2) / ((y0 - y1) + stII_c), relation) else 0,
    "ir_utils"
)
########## END - Y for calculating relations ###################


########## START - Z for counting As ###################
ir.def_constraint_class(
    'IncrementZIfA',
    lambda x, pos, var: [*var([('X',x)]), *var([('Z',pos-1), ('Z',pos)])],
    lambda nt, z0, z1: 1 if (nt == 0 and z0 + 1 == z1) or (nt != 0 and z0 == z1) else 0,
    "ir_utils"
)
ir.def_constraint_class(
    'ZSmallerThan',
    lambda i, A, var: var([('Z',i-1)]),
    lambda z, A: 1 if z <= A else 0,# p1 if x <= cutOff else 0 print(x, c)
    "ir_utils"
)
########## End - Z for counting As ###################


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
    if range_[0] <= round(x, 3) <= range_[1]:
        return 1
    return 0


def bpenergy(x, y, is_terminal=False):
    if x == 4 or y == 4:
        return 0
    bpidx = bpindex_tab[x][y]
    return (params_bp[bpidx//2 + (3 if is_terminal else 0)]
            if bpidx >= 0 else -math.inf)


def hairpin_position_finding(loop, seq):
    c = len(seq[loop[0]:loop[1]]) -seq[loop[0]:loop[1]].count('X')
    p1 = seq[:loop[0]].count('X')
    p2 = seq[:loop[1]].count('X')
    return [p2, p1], c


def stem_position_finding(stem, seq):
    p = []
    for s in stem:
        c = len(seq[s[0]:s[1]]) -seq[s[0]:s[1]].count('X')
        p.append(seq[:s[1]].count('X'))
        p.append(seq[:s[0]].count('X'))
    return p, c

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

#min max
relations_range = {
        'st2_hl2': [0.6, 3.4],
        'st3_hl3': [0.6, 3.5],
        'sI_sII': [1.6, 2.9],
        'sII_sIII': [0.5, 0.8],
        'sI_sIII': [1.2, 1.5]
    }



def create_model(model_input, feature_weights):

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


    
    # ss_alpha = model_input.structure_span['gaps']['alpha']
    ss_beta = model_input.structure_span['gaps']['beta']
    ss_gamma = model_input.structure_span['gaps']['gamma']

    # min_alpha = 2
    # max_alpha = 8
    min_beta = 1
    max_beta = 6
    min_gamma = 6
    max_gamma = 8

    # # cur_min_alpha = 22
    # # limit alpha to be between cur_min_alpha + 2 and cur_min_alpha + 8 --> 24 and 30
    # pos_A = possible_gaps.index(ss_alpha[0]) + 1     # plus one because we have an artificial zero in variable Y
    # pos_B = possible_gaps.index(ss_alpha[1]) + 1     # plus one because we have an artificial zero in variable Y
    # model.add_constraints(DiffInRange(pos_A, pos_B, min_alpha, max_alpha, var, 'Y'))

    # cur_min_beta = 9
    # limit beta to be between cur_min_beta + 1 and cur_min_beta + 6 --> 10 and 15
    pos_A = possible_gaps.index(ss_beta[0])      # plus one because we have an artificial zero in variable Y
    pos_B = possible_gaps.index(ss_beta[1]) + 1     # plus one because we have an artificial zero in variable Y
    model.add_constraints(DiffInRange(pos_A, pos_B, min_beta, max_beta, var, 'Y'))

    # cur_min_gamma = 12
    # limit gamma to be between cur_min_gamma + 6 and cur_min_gamma + 8 --> 18 and 20
    pos_A = possible_gaps.index(ss_gamma[0])      # plus one because we have an artificial zero in variable Y
    pos_B = possible_gaps.index(ss_gamma[1]) + 1    # plus one because we have an artificial zero in variable Y
    model.add_constraints(DiffInRange(pos_A, pos_B, min_gamma, max_gamma, var, 'Y'))

    # ## constraints for relation
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']
    
    stem_beta = stem_position_finding(stems[1], iupac)
    stem_gamma = stem_position_finding(stems[2], iupac)

    hl_beta = hairpin_position_finding(loops['hl2'], iupac)
    hl_gamma = hairpin_position_finding(loops['hl3'], iupac)

    model.add_constraints(CalculateRelationStemHairpin(stem_beta[0], stem_beta[1], hl_beta[0], hl_beta[1], 'st2_hl2', var))
    model.add_constraints(CalculateRelationStemHairpin(stem_beta[0], stem_beta[1], hl_beta[0], hl_beta[1], 'st2_hl2', var))
    model.add_constraints(CalculateRelationStemHairpin(stem_gamma[0], stem_gamma[1], hl_gamma[0], hl_gamma[1], 'st3_hl3', var))

    sII = hairpin_position_finding([stems[1][0][0], stems[1][1][1]], iupac)
    sIII = hairpin_position_finding([stems[2][0][0], stems[2][1][1]], iupac)

    # stem I
    a_,b_ = stem_position_finding(stems[0][0], iupac)
    c_,d_  = stem_position_finding(stems[0][1], iupac)
    sI_c = b_ + d_
    sI = [c_[2], c_[1],  a_[2]], sI_c

    model.add_constraints(CalculateRelationStructure(sII[0], sII[1], sIII[0], sIII[1], 'sII_sIII', var))
    model.add_constraints(CalculateRelationStructureI(sI[0], sI[1], sII[0], sII[1], 'sI_sII', var))
    model.add_constraints(CalculateRelationStructureI(sI[0], sI[1], sIII[0], sIII[1], 'sI_sIII', var))



    if target_length:
        n_no_gaps = (len(possible_gaps) + target_length) - n
        model.add_constraints(CheckLastY(n_y, n_no_gaps, var))
    else:
        # Total Length
        model.add_functions([NotGap(i) for i in range(n)], 'totLength')
        model.set_feature_weight(feature_weights['totLength'].weight, 'totLength')


    # ######## additional constraints from simulation analysis ########
    # loops = model_input.structure_span['loop']
    # upk1 = range(*loops['upk1'])
    # n_z = len(upk1) + 1
    # model.add_variables(n_z, n_z + 1, name='Z')

    # model.restrict_domains([('Z',0)], (0, 0))
    # for i in range(1, n_z):
    #     model.restrict_domains([('Z',i)], (0, i-1))
    # var = model.idx

    # for i, pos in enumerate(upk1, start=1):
    #     model.add_constraints(IncrementZIfA(pos, i, var))
    # max_A = 8
    # model.add_constraints(ZSmallerThan(n_z, max_A, var))

    return model

def create_model_target(model_input):

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

    # ss_alpha = model_input.structure_span['gaps']['alpha']
    ss_beta = model_input.structure_span['gaps']['beta']
    ss_gamma = model_input.structure_span['gaps']['gamma']

    # min_alpha = 2
    # max_alpha = 8
    min_beta = 1
    max_beta = 6
    min_gamma = 6
    max_gamma = 8

    # cur_min_beta = 9
    # limit beta to be between cur_min_beta + 1 and cur_min_beta + 6 --> 10 and 15
    pos_A = possible_gaps.index(ss_beta[0])     # plus one because we have an artificial zero in variable Y
    pos_B = possible_gaps.index(ss_beta[1]) + 1     # plus one because we have an artificial zero in variable Y
    model.add_constraints(DiffInRange(pos_A, pos_B, min_beta, max_beta, var, 'Y'))

    # cur_min_gamma = 12
    # limit gamma to be between cur_min_gamma + 6 and cur_min_gamma + 8 --> 18 and 20
    pos_A = possible_gaps.index(ss_gamma[0])     # plus one because we have an artificial zero in variable Y
    pos_B = possible_gaps.index(ss_gamma[1]) + 1     # plus one because we have an artificial zero in variable Y
    model.add_constraints(DiffInRange(pos_A, pos_B, min_gamma, max_gamma, var, 'Y'))

    # ## constraints for relation
    stems = model_input.structure_span['stem']
    loops = model_input.structure_span['loop']
    
    stem_beta = stem_position_finding(stems[1], iupac)
    stem_gamma = stem_position_finding(stems[2], iupac)

    hl_beta = hairpin_position_finding(loops['hl2'], iupac)
    hl_gamma = hairpin_position_finding(loops['hl3'], iupac)


    model.add_constraints(CalculateRelationStemHairpin(stem_beta[0], stem_beta[1], hl_beta[0], hl_beta[1], 'st2_hl2', var))
    model.add_constraints(CalculateRelationStemHairpin(stem_beta[0], stem_beta[1], hl_beta[0], hl_beta[1], 'st2_hl2', var))
    model.add_constraints(CalculateRelationStemHairpin(stem_gamma[0], stem_gamma[1], hl_gamma[0], hl_gamma[1], 'st3_hl3', var))

    sII = hairpin_position_finding([stems[1][0][0], stems[1][1][1]], iupac)
    sIII = hairpin_position_finding([stems[2][0][0], stems[2][1][1]], iupac)

    # stem I
    a_,b_ = stem_position_finding(stems[0][0], iupac)
    c_,d_  = stem_position_finding(stems[0][1], iupac)
    sI_c = b_ + d_
    sI = [c_[2], c_[1],  a_[2]], sI_c

    model.add_constraints(CalculateRelationStructure(sII[0], sII[1], sIII[0], sIII[1], 'sII_sIII', var))
    model.add_constraints(CalculateRelationStructureI(sI[0], sI[1], sII[0], sII[1], 'sI_sII', var))
    model.add_constraints(CalculateRelationStructureI(sI[0], sI[1], sIII[0], sIII[1], 'sI_sIII', var))

    if target_length:
        n_no_gaps = (len(possible_gaps) + target_length) - n
        model.add_constraints(CheckLastY(n_y, n_no_gaps, var))
    else:
        # Total Length
        model.add_functions([NotGap(i) for i in range(n)], 'totLength')

    return model

