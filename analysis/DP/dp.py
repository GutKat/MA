#!/usr/bin/env python

# --------------------------------------------------------------------
# dp.py
#
# Deformation Profile's main script
#
# history:
#    20090825 - 1.0.0 - JAC - first version
# --------------------------------------------------------------------

import os
import sys

import dp_2d
import dp_lib
import dp_match
import dp_util

E_ = enumerate
X_ = lambda l: range( len(l) )


class Command:
    def __init__(self):
        self.out_dir = ""
        
        self.ref_pdb = ("", 0)
        self.cmp_pdbs = []
        
        self.alignments = []
        self.squares = []
        
        self.save_matrix = True
        self.save_svg = True
        
    def run(self):
        self.parse_input()
        self.check_input()

        if( not (self.save_matrix or self.save_svg) ):
            dp_util.Msg.fatal( "Nothing to do!\nCheck 'save_matrix' and 'save_svg' parameters." )
        
        # starts the match class
        match = dp_match.Match( self.alignments )
        
        # opens reference file
        dp_util.Msg.out( "opening reference file: '%s'\n" %self.ref_pdb[0] )
        match.set_reference( self.ref_pdb )
        
        dp = dp_lib.DeformationProfile( match, self.squares )
        
        for cmp_pdb in self.cmp_pdbs:
            # opens comparing model
            dp_util.Msg.out( "opening comparing file: '%s'\n" %cmp_pdb[0] )
            set_comparing = match.set_comparing( cmp_pdb )
            
            # computes 'deformation profile'
            dp_util.Msg.out( "comparing models...\n" )
            
            match.show( os.path.basename(self.ref_pdb[0]), os.path.basename(cmp_pdb[0]) )
            dp.compute( )
            
            # saves a data file
            if( self.save_matrix ):
                dp_util.Msg.out( "saving data file...\n" )
                dp.matrix_save( self.file_base_name( cmp_pdb[0] ) + ".dat" )
                 
            # saves an SVG file
            if( self.save_svg ):
                dp_util.Msg.out( "saving svg file...\n" )
                dp.svg_save( self.file_base_name( cmp_pdb[0] ) + ".svg" )
            
    def file_base_name(self, pdb_name ):
        if( self.out_dir != "" ):
            pdb_name = "%s%s%s" %(self.out_dir, os.sep, os.path.basename(pdb_name))
            
        if( pdb_name.endswith( ".pdb" ) ):
            bname = pdb_name[:-4]
        else:
            bname = pdb_name

        return( bname )

    def parse_input( self ):
        if( len(sys.argv) != 3 ):
            dp_util.Msg.usage()

        if( sys.argv[1] == "-c" ):
            cfg = sys.argv[2]

            if( os.path.isfile( cfg ) ):
                # initializes all config variables
                config_namespace = {
                    'out_dir': '.',
                    'ref_model': '',
                    'cmp_model': '',
                    'cmp_list': '',
                    'aligns': [],
                    'helices': [],
                    'loops': [],
                    'draw': [],
                    'matrix': True,
                    'svg' : True,
                    'quiet_err': False,
                    'quiet_out': False,
                }

                # calls config file
                config = open( cfg ).read()
                exec(config, globals(), config_namespace)

                # translates data
                self.out_dir = config_namespace['out_dir']
                
                self.ref_pdb = config_namespace['ref_model']

                if( config_namespace['cmp_model'] != "" ):
                    self.cmp_pdbs = config_namespace['cmp_model']
                
                if( config_namespace['cmp_list'] != "" ):
                    self.cmp_pdbs = []
                    self.parse_cmp_list( config_namespace['cmp_list'] )

                self.alignment = []
                for (ref_chain, ref_start, cmp_chain, cmp_start, length) in config_namespace['aligns'] :
                    self.alignments.append( dp_match.Alignment(length, ref_chain, ref_start, cmp_chain, cmp_start) )
                    
                dp_util.Msg.STDERR_QUIET = config_namespace['quiet_err']
                dp_util.Msg.STDOUT_QUIET = config_namespace['quiet_out']
                
                # secondary structure definition
                index_aux = {}
                
                hs = []
                for (name, first5, length5, first3, length3) in config_namespace['helices']:
                    index_aux[name] = "H"
                    hs.append( dp_2d.Helix( name, first5, length5, first3, length3 ) )
                
                ls = []
                for (name, first, length) in config_namespace['loops'] :
                    index_aux[name] = "L"
                    ls.append( dp_2d.Loop( name, first, length ) )
                
                ss = dp_2d.SecondaryStructure( hs, ls )
                
                for key in config_namespace['draw']:
                    self.squares.extend( self.parse_draw_key( ss, key, index_aux ) )
                
                # action
                self.save_matrix = config_namespace['matrix']
                self.save_svg = config_namespace['svg']
            else:
                dp_util.Msg.fatal( "'%s' file not found\n" %cfg )
        elif( sys.argv[1] == "-o" ):
            pdb = sys.argv[2]
            
            if( os.path.isfile( pdb ) ):
                dp_util.show_data( pdb )
                quit()
            else:
                dp_util.Msg.fatal( "'%s' file not found\n" %pdb )
        else:
            # minimal usage mode
            self.ref_pdb = (sys.argv[1], 0)
            self.cmp_pdbs = [(sys.argv[2], 0)]

    # performs some obvious check in the input
    def check_input(self):
        print(self.cmp_pdbs)
        if( (self.out_dir) != "" and (not os.path.isdir( self.out_dir )) ):
            dp_util.Msg.fatal( "Output directory not found: '%s'\ncheck 'out_dir' parameter" %self.out_dir )

        if( not os.path.isfile( self.ref_pdb[0] ) ):
            dp_util.Msg.fatal( "Reference file not found: '%s'\ncheck 'ref_model' parameter" %self.ref_pdb[0] )

        for (fname, model) in self.cmp_pdbs:
            if( not os.path.isfile( fname ) ):
                dp_util.Msg.fatal( "Comparing file not found: '%s'\ncheck 'cmp_model' or 'cmp_list'  parameter" %fname )

    def parse_cmp_list(self, cmp_list):
        if( os.path.isfile(cmp_list) ):
            fi = open( cmp_list, "r" )
            
            for line in fi:
                line = line.strip()
                if( (line == "") or (line.startswith( "#" )) ):
                    continue
                
                data = line.split(";")
                if( len(data) == 2 ):
                    if( data[1].isdigit() ):
                        self.cmp_pdbs.append( (data[0], int(data[1])) )
                    else:
                        dp_util.Msg.fatal( "Model number must be a number in line: %s" %line)
                else:
                    dp_util.Msg.fatal( "Expecting: 'file name;# model' in line: %s" %line)
        else:
            dp_util.Msg.fatal( "Comparing list file not found: '%s'\ncheck 'cmp_list' parameter" %cmp_list )

    def parse_draw_key(self, ss, key, index_aux):
        result = []
        strip = lambda x: x.strip()
        
        name = ""
        if( ":" in key ):
            (data, name) = map( strip, key.split( ":" ) )
        else:
            data = key

        if( "x" in data ):
            (s1, s2) = map( strip, data.split( "x" ) )
            if( index_aux.has_key(s1) and index_aux.has_key(s2) ):
                sqr_type = index_aux[s1] + index_aux[s2]
                if( sqr_type == "HH" ):
                    result = ss.square_hh( s1, s2, name )
                elif( sqr_type == "LL" ):
                    result = ss.square_ll( s1, s2, name )
                elif( sqr_type == "HL" ):
                    result = ss.square_hl( s1, s2, name )
                elif( sqr_type == "LH" ):
                    result = ss.square_lh( s1, s2, name )
            else:
                dp_util.Msg.fatal( "Syntax error in draw key: '%s'\ncheck 'draw' parameter" %key )
        else:
            s = key.strip()
            
            if( index_aux.has_key(s) ):
                if( index_aux[s] == "H" ):
                    result = ss.square_helix( s, name )
                elif( index_aux[s] == "L" ):
                    result = ss.square_loop( s, name )
            else:
                dp_util.Msg.fatal( "Syntax error in draw key: '%s'\ncheck 'draw' parameter" %key )
        
        return( result )


#
# MAIN
#
cmd = Command()
cmd.run()

# import subprocess
# subprocess.run(['convert', '/Users/katringutenbrunner/Desktop/MA/working/opt/Deformation_Profile/examples/ex1/b.svg', '/Users/katringutenbrunner/Desktop/MA/working/opt/Deformation_Profile/examples/ex1/b.png'])