# --------------------------------------------------------------------
# dp.py
#
# Miscelaneous utility classes
#
# history:
#    20090825 - 1.0.0 - JAC - first version
# --------------------------------------------------------------------

import os
import sys

from Bio.PDB import *

class Msg:
    STDERR_QUIET = False
    STDOUT_QUIET = False
    
    def fatal( msg ):
        sys.stderr.write( "Fatal Error!\n%s\n" %msg )
        quit()
        
    def out( msg, stderr=True ):
        if( stderr ):
            if( not Msg.STDERR_QUIET ):
                sys.stderr.write( msg )
                sys.stderr.flush()
        else:
            if( not Msg.STDOUT_QUIET ):
                sys.stdout.write( msg )
        
    def usage():
        cmd = os.path.basename(sys.argv[0])
        
        sys.stderr.write( "\n-------------------\n" )
        sys.stderr.write( "Deformation Profile\n" )
        sys.stderr.write( "-------------------\n" )
        sys.stderr.write( "\nUsage:\n" )
        sys.stderr.write( "\t%s <reference pdb> <comparing pdb>\n" %cmd )
        sys.stderr.write( "\t%s -c <config_file>\n" %cmd )
        sys.stderr.write( "\t%s -o <pdb_file>\n\n" %cmd )
        sys.stderr.write( "For more information see 'dps_manual.pdf'\n\n" )
        quit()
        
    fatal = staticmethod( fatal )
    out = staticmethod( out )
    usage = staticmethod( usage )

#
# This implements a short version of S-W local alignment algorithm to obtain the longest gap-free sequence.
#
class SeqAligner:
    def __init__(self, match=2, miss=-1,indel=-4):
        self._match = match
        self._miss = miss
        self._indel = indel
        
        self.length = 0
        self.start1 = 0
        self.start2 = 0
        self.bases1 = []
        self.bases2 = []
        
    def matrix_build( self, s1, s2 ):
        F = [[0] * (len(s2)+1) for i in range(len(s1)+1)]
      
        max_coords = (0, 0)
        for i in range(1, len(s1)+1):
            for j in range(1, len(s2)+1):
                if( s1[i-1] == s2[j-1] ):
                    score = self._match
                else:
                    score = self._miss
                     
                choice0 = 0
                choice1 = F[i-1][j-1] + score
                choice2 = F[i-1][j] + self._indel
                choice3 = F[i][j-1] + self._indel
                
                F[i][j] = max( choice1, choice2, choice3 )
                
                if( F[max_coords[0]][max_coords[1]] < F[i][j] ):
                    max_coords = (i, j)
        
        return (F, max_coords)
    
    def align( self, s1, s2 ):
        self.length = 0
        self.start1 = 0
        self.start2 = 0
        self.bases1 = []
        self.bases2 = []
        
        (F, (mi, mj)) = self.matrix_build( s1, s2 )
        
        i = mi
        j = mj
            
        while (i > 0 and j > 0):
            if( (F[i][j] >= F[i-1][j]) and (F[i][j] >= F[i][j-1]) ):
                i -= 1
                j -= 1
            else:
                break
        
        self.length = mi-i
        self.start1 = i
        self.start2 = j
        self.bases1 = s1[i:mi]
        self.bases2 = s2[j:mj]    
            
            
#
# PDB utility functions
#
def get_sequence( chain ):
    return( [r.get_resname().strip() for r in chain] )

def get_model( pdb, model_num ):
    parser = PDBParser()
    
    struct = parser.get_structure( "X", pdb )
    
    if( model_num < len( struct ) ):
        return( struct[model_num] )
    else:
        Msg.out( "No model '%d' in pdb file '%s'" %(model_num, pdb) )
        
def show_data(pdb):
    parser = PDBParser()
    
    struct = parser.get_structure( "X", pdb )
    
    Msg.out( "PDB '%s' (%d models):\n" %(pdb, len(struct)), False )
    
    for model in struct:
        Msg.out( "\tModel '%s' (%d chains):\n" %(model.get_id(), len(model) ), False )
        
        for chain in model:
            Msg.out( "\t\tChain '%s' (%s residues):\n" %(chain.get_id(), len(chain)), False )
            
            txt1, txt2 = "", ""
            for residue in chain:
                txt1 += "%4d " %residue.get_id()[1]
                txt2 += "%4s " %residue.get_resname()
            
            Msg.out( "\t\t%s\n\t\t%s\n\n" %(txt1, txt2), False )

if __name__ == '__main__':
    s1 =     "GCUXAGAUXAYYYGUXUAXXXX"
    s2 = "AUUGGCUUAGAUCAAGUGUAGUAUCUGUUCUUUUCAGUGUA"

    sa = SeqAligner()
    sa.align( s1, s2 )
    
    print(sa.length)
    print(sa.start1)
    print(sa.start2)
    print(sa.bases1)
    print(sa.bases2)
