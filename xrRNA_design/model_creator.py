import infrared as ir
import RNA
import infrared.rna as rna
from collections import namedtuple
import utils as ut

ModelInput = namedtuple("ModelInput", "structures anti_structures iupac var_region_simple var_region_mixed")


def create_model(model_input):


    n = len(model_input.target_length)

    model = ir.Model()


    #First n variables descibe the base positions. Have 4 different bases and a gap (domain of five)
    model.add_variables(n,5)

    #Variables for the variable region.
    model.add_variables(n,n)

    for target in model_input.targets:
        ss = rna.parse(target)
        model.add_constraints(rna.BPComp(i, j) for (i, j) in ss)
        model.add_functions([rna.BPEnergy(i, j, (i-1, j+1) not in ss)])

    
    for regions in model_input.var_stem_regions:
        pass

        
    for start,stop in model_input.var_loop_regions:
        model.add_constraints(EmptyIfLeftEmpty(i) for i in range(var_simple[0],var_simple[1]))

    return model


