import sys,os
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

from RNA_assessment import RNA_normalizer
from operator import attrgetter
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()

def create_index(file, start=1, length=None):
	structure = parser.get_structure("1", file)
	model = structure[0]
	chain = model['A']
	file_index = file.replace('.pdb', '.index')
	if not length:
		length = len(chain) - start
	output_index = f'A:{start}:{length}'
	with open(file_index, 'w') as f:
		f.write(output_index)

	return file_index

RESIDUES_LIST = "/scr/aldea/kgutenbrunner/working/analysis/RNA_assessment/data/residues.list"
ATOMS_LIST = "/scr/aldea/kgutenbrunner/working/analysis/RNA_assessment/data/atoms.list"

def CleanFormat(f):
	os.system( "mac2unix -q %s" %f )
	os.system( "dos2unix -q %s" %f )
	

def normalize_structure(struct, out_file = None, index_file=None, extract_file = None):
	pdb_normalizer = RNA_normalizer.PDBNormalizer( RESIDUES_LIST, ATOMS_LIST )
	ok = pdb_normalizer.parse( struct, out_file )
	if not ok:
		sys.stderr.write("ERROR: structure not normalized!\n")
	if not extract_file is None:
		coords=open(index_file).read()
		RNA_normalizer.extract_PDB(SOLUTION_NORMAL,coords, extract_file)
		sys.stderr.write("INFO:	structure extracted\n")
		


# PVALUE set according to Hajdin et al., RNA (7) 16, 2010, either "+" or "-"
def calc_RMSD(res_struct, sol_struct, PVALUE = "-"):	

	sol_raw_seq = sol_struct.raw_sequence()
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	# sys.stderr.write("INFO Partial RMSD --> %f\n" %rmsd )
	pvalue = comparer.pvalue( rmsd, len(sol_raw_seq), PVALUE )
	# sys.stderr.write("INFO Partial P-Value --> %e\n" %pvalue )
	return(rmsd, pvalue)


def InteractionNetworkFidelity(res_struct, sol_struct):
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	INF_ALL = comparer.INF( sol_struct, res_struct, type="ALL" )
	DI_ALL = rmsd / INF_ALL
	INF_WC = comparer.INF( sol_struct, res_struct, type="PAIR_2D" )
	INF_NWC = comparer.INF( sol_struct, res_struct, type="PAIR_3D" )
	INF_STACK = comparer.INF( sol_struct, res_struct, type="STACK" )
	return (rmsd,DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK)
