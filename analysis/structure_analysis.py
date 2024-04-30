import os
import argparse
import subprocess
import utils as ut
import re
import pandas as pd
import barnaba as bb
from pathlib import Path
from tqdm import tqdm


from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

parser=argparse.ArgumentParser()
# Define mutually exclusive group
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-t', '--targets', nargs='+', help='List of target files')
group.add_argument('-f', '--folder', help='directory containing target files, also used for all output files')

# get all the arguments of the command

parser.add_argument('-r', '--reference', help='reference file', required=True)
parser.add_argument('-o', '--output', help='Output folder for storing output and all intermediate created files')
# parser.add_argument('-s', '--scores', help='scores which should be calculated')
args=parser.parse_args()

# Check if targets are provided and require output folder
if args.targets and not args.output:
    parser.error("Output folder is required when specifying targets.")
    raise TypeError

# create complete paths from the arguments
pwd = os.getcwd()
reference = args.reference #pwd + '/' +


if args.folder:
    input_folder = args.folder
    # if pwd not in input_folder:
      #   input_folder = pwd + '/' + input_folder
        ### not good --> need to change this
    targets = [input_folder + '/' + file for file in os.listdir(input_folder) if file.endswith('pdb')]
    output_folder = input_folder
else:
    targets = args.targets
    output_folder = args.output
    # if pwd not in output_folder:
    #    output_folder = pwd + '/' + output_folder

output_file = output_folder + '/' + 'result.csv'

# print metadata from input to user
print('-' * 100)
print('Input data')
metadata_text =  f'Reference file: {os.path.basename(reference)}\nTarget files: {[os.path.basename(target) for target in targets]}\nOutput files are stored in: {output_folder}\nResults stored in: {output_file}'
print(metadata_text)
print('-' * 100)



columns = ['Name', 'RMSD', 'p-value', 'DI_ALL', 'INF_ALL', 'INF_WC', 'INF_NWC','INF_STACK', 'CAD-Score', 'TM-score', 'GTD TS', 'eRMSD']
df = pd.DataFrame(columns=columns)


norm_reference = reference.replace('.pdb', '_norm.pdb')
ut.normalize_structure(reference, out_file=norm_reference)
reference_index = ut.create_index(norm_reference, 3, 87)

voronota = '/Users/katringutenbrunner/Desktop/MA/working/opt/rna_analysis/voronota/voronota-cadscore'
# voronota = '/scr/aldea/kgutenbrunner/opt/rna_analysis/voronota_1.28.4132/voronota-cadscore'


TM_regexr = r'TM-score *= *(\d+.\d+)'
GTD_TS_regexr = r'GDT-TS-score *= *(\d+.\d+)'

TM = '/Users/katringutenbrunner/Desktop/MA/working/opt/rna_analysis/TMscore'
#TM = '/scr/aldea/kgutenbrunner/opt/rna_analysis/TMscore'

print('Starting metric-calculations')
for target in (targets): #tqdm
    target_name = Path(target).stem
    norm_target = target.replace('.pdb', '_norm.pdb')
    ut.normalize_structure(target, out_file=norm_target)
    targets_index = ut.create_index(target, 4, 87)
    try:
    # calculate RMSD for RNA structures
        rmsd, pval = ut.calc_RMSD(native_file=norm_reference, native_index=reference_index, prediction_file=norm_target, prediction_index=targets_index)
    # calculate InteractionNetworkFidelity and Deformation Index for RNA structures
        rmsd2, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK = ut.InteractionNetworkFidelity(native_file=norm_reference, native_index=reference_index, prediction_file=norm_target, prediction_index=targets_index)
    # print('RMSD: %.3f; P-value: %.3e; Deformation Index: %.3f; INF_all: %.3f; INF_wc: %.3f; INF_nwc: %.3f; INF_stack: %.3f'%(rmsd,pval,DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK))
    except TypeError:
        print("RMSD, DI, INF could not be calculated because sequences are not of the same length")
        rmsd = pval=DI_ALL= INF_ALL = INF_WC=INF_NWC=INF_STACK = None


    
    cmd_vor = f'{voronota} -m {reference} -t {target}'
    voronota_result = subprocess.check_output(cmd_vor, shell=True, text=True)
    voronota_result = voronota_result.replace('\n','').split(' ')
    cad_score = float(voronota_result[4])
    area_target = float(voronota_result[5])
    area_model = float(voronota_result[6])


    cmd_TM = f'{TM} {reference} {target}'
    TM_result = subprocess.check_output(cmd_TM, shell=True, text=True)
    try:
        TM_score = float(re.search(TM_regexr, TM_result)[1])
        GTD_TS = float(re.search(GTD_TS_regexr, TM_result)[1])
    except:
        print("TM-Score did not work. TM score and GTD_TS were not calculated")
        TM_score = None
        GTD_TS = None
    try:
        ermsd = bb.ermsd(reference,target)[0]
    except AssertionError:
        print("ermsd could not be calculated because sequences are not of the same length")
        ermsd = None

    results = [target_name, rmsd, pval, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK, cad_score, TM_score, GTD_TS, ermsd]
    df.loc[len(df.index)] = results
print('Finished metric-calculations')

print('-' * 100)

print(df)
df.to_csv(output_file)
print(f'\nResults were stored in {output_file}\n')
