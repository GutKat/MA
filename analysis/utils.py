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

RESIDUES_LIST = "/Users/katringutenbrunner/Desktop/MA/working/analysis/RNA_assessment/data/residues.list"
ATOMS_LIST = "/Users/katringutenbrunner/Desktop/MA/working/analysis/RNA_assessment/data/atoms.list"
# RESIDUES_LIST = "/scr/aldea/kgutenbrunner/working/analysis/RNA_assessment/data/residues.list"
# ATOMS_LIST = "/scr/aldea/kgutenbrunner/working/analysis/RNA_assessment/data/atoms.list"

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
def calc_RMSD(native_file, native_index, prediction_file, prediction_index, PVALUE = "-"):
	res_struct = RNA_normalizer.PDBStruct()
	res_struct.load( native_file, native_index )
	res_raw_seq = res_struct.raw_sequence()
	
	sol_struct = RNA_normalizer.PDBStruct()
	sol_struct.load( prediction_file, prediction_index )
	sol_raw_seq = sol_struct.raw_sequence()
	
	if( sol_raw_seq != res_raw_seq ):
		sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
		sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
		sys.stderr.write("DATA Result sequence   --> '%s'\n" %res_raw_seq )
		return(-1)
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	# sys.stderr.write("INFO Partial RMSD --> %f\n" %rmsd )
	pvalue = comparer.pvalue( rmsd, len(sol_raw_seq), PVALUE )
	# sys.stderr.write("INFO Partial P-Value --> %e\n" %pvalue )
	return(rmsd, pvalue)


def InteractionNetworkFidelity(native_file, native_index, prediction_file, prediction_index):
	res_struct = RNA_normalizer.PDBStruct()
	res_struct.load( native_file, native_index )
	res_raw_seq = res_struct.raw_sequence()
	
	sol_struct = RNA_normalizer.PDBStruct()
	sol_struct.load( prediction_file, prediction_index )
	sol_raw_seq = sol_struct.raw_sequence()
	
	if( sol_raw_seq != res_raw_seq ):
		sys.stderr.write("ERROR Result sequence != Solution sequence!\n")
		sys.stderr.write("DATA Solution sequence --> '%s'\n" %sol_raw_seq )
		sys.stderr.write("DATA Result sequence   --> '%s'\n" %res_raw_seq )
		# return(-1)
	# computes the RMSD
	comparer = RNA_normalizer.PDBComparer()
	rmsd = comparer.rmsd( sol_struct, res_struct )
	INF_ALL = comparer.INF( sol_struct, res_struct, type="ALL" )
	DI_ALL = rmsd / INF_ALL
	INF_WC = comparer.INF( sol_struct, res_struct, type="PAIR_2D" )
	INF_NWC = comparer.INF( sol_struct, res_struct, type="PAIR_3D" )
	INF_STACK = comparer.INF( sol_struct, res_struct, type="STACK" )
	return (rmsd,DI_ALL, INF_ALL, INF_WC, INF_NWC,INF_STACK)