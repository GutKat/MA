import os
import argparse
import subprocess
import utils as ut
import re
import pandas as pd
import barnaba as bb
from pathlib import Path
from RNA_assessment import RNA_normalizer


from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

parser=argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-t', '--targets', nargs='+', help='List of target files')
group.add_argument('-f', '--folder', help='directory containing target files, also used for all output files')


# parser.add_argument('-r', '--reference', help='reference file', required=True)
parser.add_argument('-o', '--output', help='Output folder for storing output and all intermediate created files')
args=parser.parse_args()

if args.targets and not args.output:
    parser.error("Output folder is required when specifying targets.")


# create complete paths from the arguments
pwd = os.getcwd()
reference =  '/scr/aldea/kgutenbrunner/data/tbfv_sim/DTV_xrRNA1_pk1_new_pk3_new_seq/cluster1_min.pdb'
truncated_ref = '/scr/aldea/kgutenbrunner/data/tbfv_sim/DTV_xrRNA1_pk1_new_pk3_new_seq/edits/truncated.pdb'


if args.folder:
    input_folder = args.folder
    targets = [input_folder + '/' + file for file in os.listdir(input_folder) if file.endswith('pdb')]
    output_folder = input_folder
else:
    targets = args.targets
    output_folder = args.output

output_file = output_folder + '/' + 'result.csv'


columns = ['Name', 'RMSD', 'p-value', 'DI_ALL', 'INF_ALL', 'INF_WC', 'INF_NWC','INF_STACK', 'CAD-Score', 'TM-score to ref','TM-score to pred', 'eRMSD']
df = pd.DataFrame(columns=columns)

norm_reference = reference.replace('.pdb', '_norm.pdb')
ut.normalize_structure(reference, out_file=norm_reference)
reference_index = ut.create_index(norm_reference, 3, 87)
res_struct = RNA_normalizer.PDBStruct()
res_struct.load( norm_reference, reference_index )


deleting = [norm_reference, reference_index]

TM_regexr = r'TM-score= (\d+.\d+)'
TM = '/scr/aldea/kgutenbrunner/opt/rna_analysis/RNAalign/RNAalign'
TM_alignment = ['>reference', 'CAGGGGGUGAUGUGGCAGCGCACCACGACAUCGUGACGGGAAGAGGUCGUCCCCGACGCAUCAUCUCUCUAGGGCAUUUUCGUGAGACCCUCAU', '>comparison', '']
MCQ = '/scr/aldea/kgutenbrunner/opt/rna_analysis/mcq4structures/mcq-cli/mcq-local'
VORONOTA = '/scr/aldea/kgutenbrunner/opt/rna_analysis/voronota_1.28.4132/voronota-cadscore'

mcq_folder = output_folder + '/mcq'
os.mkdir(mcq_folder)

name_pattern = r'(\d+_\d+)\/.*\/(clust\d+)'

def main():
    # print metadata from input to user
    print('-' * 100)
    print('Input data')
    metadata_text =  f'Reference file: {os.path.basename(reference)}\nTarget files: {[os.path.basename(target) for target in targets]}\nOutput files are stored in: {output_folder}\nResults stored in: {output_file}'
    print(metadata_text)
    print('-' * 100)

    
    print('Starting metric-calculations')
    for target in (targets): #tqdm
        target_name = re.search(name_pattern,target)[1] + '_' + re.search(name_pattern,target)[2]


        
        norm_target = target.replace('.pdb', '_norm.pdb')
        ut.normalize_structure(target, out_file=norm_target)
        target_index = ut.create_index(target, 1, 87)

        sol_struct = RNA_normalizer.PDBStruct()
        sol_struct.load( norm_target, target_index )
        sol_seq = sol_struct.raw_sequence()

        try:
        # calculate RMSD for RNA structures
            rmsd, pval = ut.calc_RMSD(res_struct=res_struct,  sol_struct=sol_struct)
        # calculate InteractionNetworkFidelity and Deformation Index for RNA structures
            rmsd2, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK = ut.InteractionNetworkFidelity(res_struct=res_struct, sol_struct=sol_struct)
        except TypeError:
            print("RMSD, DI, INF could not be calculated because sequences are not of the same length")
            rmsd = pval=DI_ALL= INF_ALL = INF_WC=INF_NWC=INF_STACK = None


        
        cmd_vor = f'{VORONOTA} -m {reference} -t {target}'
        voronota_result = subprocess.check_output(cmd_vor, shell=True, text=True)
        voronota_result = voronota_result.replace('\n','').split(' ')
        cad_score = float(voronota_result[4])
        # area_target = float(voronota_result[5])
        # area_model = float(voronota_result[6])

        TM_alignment[-1] = f'--{sol_seq}-----'
        TM_ali_name = target.replace('pdb', 'ali')
        with open(f'{TM_ali_name}', 'w') as f:
            for line in TM_alignment:
                f.write(line)
                f.write('\n')
        cmd_TM = f'{TM} {reference} {target} -I {TM_ali_name}'
        TM_result = subprocess.check_output(cmd_TM, shell=True, text=True)
        try:
            TM_score_ref = float(re.findall(TM_regexr, TM_result)[0])
            TM_score_pred = float(re.findall(TM_regexr, TM_result)[1])
        except:
            print("TM-Score did not work. TM score and GTD_TS were not calculated")
            TM_score_ref = TM_score_pred = None

        try:
            ermsd = bb.ermsd(truncated_ref,target)[0]
        except AssertionError:
            print("ermsd could not be calculated because sequences are not of the same length")
            ermsd = None

        results = [target_name, rmsd, pval, DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK, cad_score, TM_score_ref, TM_score_pred, ermsd]
        df.loc[len(df.index)] = results
        deleting.append(target_index)
        deleting.append(norm_target)
        deleting.append(f'{target}.mcout')
        deleting.append(f'{norm_target}.mcout')

    mcq_cmd = f'{MCQ} -d {mcq_folder} -t {reference} -T A:3:87 -M A:1:87 ' + ' '.join(targets)
    mcq_result = subprocess.check_output(mcq_cmd, shell=True, text=True)
    mcq_result = [line for line in mcq_result.split('\n') if line][-len(targets):]
    mcq_result = [float(x.split(' ')[1]) for x in mcq_result]
    df.insert(len(df.columns), "mcq", mcq_result, True)
    
    print('Finished metric-calculations')
    print('-' * 100)
    
    print(df)
    df.to_csv(output_file)
    print(f'\nResults were stored in {output_file}\n')

    # delete all the created files
    [os.remove(os.path.dirname(reference) + '/' + file) for file in os.listdir((os.path.dirname(reference))) if 'norm' in file or 'mcout' in file]
    [os.remove(file) for file in deleting if os.path.isfile(file)]
    print('Deleted all created intermediate-files')

# /scr/aldea/kgutenbrunner/opt/rna_analysis/mcq4structures/mcq-cli/mcq-local -d /scr/aldea/kgutenbrunner/working/analysis/outputs/test -t /scr/aldea/kgutenbrunner/working/simRNA/05/05_01/QRNAs/clust01_qrna.pdb -T A:1:45 -M A:1:45 /scr/aldea/kgutenbrunner/working/simRNA/04/04_30/QRNAs/clust01_qrna.pdb


if __name__ == "__main__":
    main()
    
