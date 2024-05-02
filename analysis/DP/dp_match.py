# --------------------------------------------------------------------
# dp_match.py
#
# These classes deal with matches between residues
#
# history:
#    20090825 - 1.0.0 - JAC - first version
# --------------------------------------------------------------------

import dp_util

E_ = enumerate
X_ = lambda l: range( len(l) )


ATOMS_FULL = ["C1'", "C1*", "C2'", "C2", "C2*", "C3'", "C3*", "C4'", "C4", "C4*", "C5'", "C5", "C5*", "C6", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2'", "O2*", "O2", "O3'", "O3*", "O4'", "O4*", "O4", "O5'", "O5*", "O6", "P", "OP1", "OP2"]
#ATOMS_FULL = ["C1'", "C1*", "C4'", "C4*", "P"]
#ATOMS_PHOSPHATE = ["P"]

NICE_NAME = lambda name: name.strip().upper().replace( "*", "'" )

class Alignment:
    def __init__(self, length, ref_chain, ref_start, cmp_chain, cmp_start):
        self.length = length
    
        self.ref_chain = ref_chain
        self.ref_start = ref_start
        
        self.cmp_chain = cmp_chain
        self.cmp_start = cmp_start

class Match:
    def __init__( self, alignments=None, atoms=None ):
        self.alignments = alignments
        
        if( atoms == None ):
            self.atoms = ATOMS_FULL
        else:
            self.atoms = atoms
        
        self.ref_model = None
        self.cmp_model = None
                
        # ref/cmp list are lists of tuples (residue, base_name, atom_list)
        self._ref_list = []
        self._cmp_list = []
        
    def show(self, ref_name, cmp_name):
        txt, rt1, rt2, rt3, ct1, ct2, ct3 = "", "", "", "", "", "", ""
        
        for i, (ref, cmp) in enumerate(zip(self._ref_list, self._cmp_list)):
            txt += "%3s " %i
            rt1 += "%3s " %ref[0].get_parent().get_id()
            rt2 += "%3s " %ref[0].get_id()[1]
            rt3 += "%3s " %ref[1]
            ct1 += "%3s " %cmp[0].get_parent().get_id()
            ct2 += "%3s " %cmp[0].get_id()[1]
            ct3 += "%3s " %cmp[1]
        
        
        dp_util.Msg.out( "%s\n" %("- " * 30), False )
        dp_util.Msg.out( "reference:    %s\n" %ref_name, False )
        dp_util.Msg.out( "comparing:    %s\n\n" %cmp_name, False )

        dp_util.Msg.out( "ref id:       %s\n" %rt2, False )
        dp_util.Msg.out( "cmp id:       %s\n\n" %ct2, False )

        dp_util.Msg.out( "ref chain:    %s\n" %rt1, False )
        dp_util.Msg.out( "cmp chain:    %s\n" %ct1, False )

        dp_util.Msg.out( "ref residues: %s\n" %rt3, False )
        dp_util.Msg.out( "cmp residues: %s\n" %ct3, False )
        dp_util.Msg.out( "align. index: %s\n" %txt, False )
        dp_util.Msg.out( "%s\n" %("- " * 30), False )

    def set_reference( self, ref_pdb ):
        # opens reference model
        self.ref_pdb = ref_pdb[0]
        self.ref_model_id = ref_pdb[1]        
        self.ref_model = dp_util.get_model( self.ref_pdb, self.ref_model_id )

    def set_comparing( self, cmp_pdb ):
        self.cmp_pdb = cmp_pdb[0]
        self.cmp_model_id = cmp_pdb[1]        
        self.cmp_model = dp_util.get_model( self.cmp_pdb, self.cmp_model_id )
        self.update()

    def get_length(self):
        return( len( self._ref_list ) )
    
    def get_residues(self, ndx):
        return self._get_generic(ndx, 0)

    def get_res_names(self, ndx):
        return self._get_generic(ndx, 1)
   
    def get_atoms(self, ndx):
        return self._get_generic(ndx, 2)
    
    def _get_generic(self, ndx, pos):
        return( self._ref_list[ndx][pos], self._cmp_list[ndx][pos] )
       
    def update( self ):
        self._ref_list = []
        self._cmp_list = []
         
        if( len(self.alignments) > 0 ):
            # if an alignment was customized
            self.update_from_alignment( )
        else:
            # infer from the models
            self.update_from_models( )
            
        # get atoms
        self.update_atoms()
            
    def update_from_models(self):
        sa = dp_util.SeqAligner()
        
        # for each chain in the reference model
        for ref_chain in self.ref_model:
            # for each chain in the comparing model
            for cmp_chain in self.cmp_model:
                # if both chains have the same name
                if( ref_chain.id == cmp_chain.id ):
                    # get both sequences
                    ref_seq = dp_util.get_sequence( ref_chain )
                    cmp_seq = dp_util.get_sequence( cmp_chain )
                    
                    # uses the best ungaped local alignment
                    sa.align( ref_seq, cmp_seq )                   
                    
                    #TODO check if both sequences have the same length!
                    for i in range( sa.length ):
                        self._ref_list.append( [ref_chain.child_list[sa.start1 + i], sa.bases1[i], None] )
                        self._cmp_list.append( [cmp_chain.child_list[sa.start2 + i], sa.bases2[i], None] )
            
    def update_from_alignment(self):
        # TODO: check length validity
        for align in self.alignments:
            if( not self.ref_model.has_id( align.ref_chain ) ):
                dp_util.Msg.fatal( "Can't find chain '%s' in reference model" %(align.ref_chain) )

            if( not self.cmp_model.has_id( align.cmp_chain ) ):
                dp_util.Msg.fatal( "Can't find chain '%s' in comparing model" %(align.cmp_chain) )
            
            for i in range( align.length ):
                pos = align.ref_start + i
                if( self.ref_model[align.ref_chain].has_id( pos ) ):
                    res = self.ref_model[align.ref_chain][pos]
                    self._ref_list.append( [res, res.get_resname().strip(), None] )
                else:
                    dp_util.Msg.fatal( "Can't find residue '%s' in reference chain '%s'" %(pos, align.ref_chain) )

                pos = align.cmp_start + i
                if( self.cmp_model[align.cmp_chain].has_id( pos ) ):
                    res = self.cmp_model[align.cmp_chain][pos]
                    self._cmp_list.append( [res, res.get_resname().strip(), None] )
                else:
                    dp_util.Msg.fatal( "Can't find residue '%s' in comparing chain '%s'" %(pos, align.cmp_chain) )


    def get_common_atoms( self, ref_residue, cmp_residue ):        
        ref_atoms = []
        cmp_atoms = []
        
        ref_atom_names = [a.name for a in ref_residue if a.name in self.atoms]
        cmp_atom_names = [a.name for a in cmp_residue if a.name in self.atoms]
        
        for ra_name in ref_atom_names:
            for ca_name in cmp_atom_names:
                if( NICE_NAME(ra_name) == NICE_NAME(ca_name) ):
                    ref_atoms.append( ref_residue[ra_name] )
                    cmp_atoms.append( cmp_residue[ca_name] )
                    break
    
        return( ref_atoms, cmp_atoms )
    
    def update_atoms(self):
        for i in X_( self._ref_list ):
            (ref_atoms, cmp_atoms) = self.get_common_atoms( self._ref_list[i][0], self._cmp_list[i][0] )
            
            self._ref_list[i][2] = ref_atoms
            self._cmp_list[i][2] = cmp_atoms