# --------------------------------------------------------------------
# dp_lib.py
#
# Deformation Profile's core classes
#
# history:
#    20090825 - 1.0.0 - JAC - first version
# --------------------------------------------------------------------

import copy
import dp_util
import numpy as np
import svg

from Bio.PDB import *

E_ = enumerate
X_ = lambda l: range( len(l) )

#
# Graphic parameters
#

SQR_SIDE = 16

# distance normalization
# *** TODO ***: these values must be set in the config file
NORMALIZE = False

# palette parameters
LIMIT_DOWN = 0.75

# *** TODO ***: these values must be set in the config file
STEPS_DOWN = 10
STEPS_UP = 10


# curve color parameters
COLOR_ROW_MEAN = "#00FF00" # green
COLOR_COL_MEAN = "#0000FF" # blue
COLOR_LOCAL_RMSD = "#FF0000" # red

class DeformationProfile:
    def __init__(self, match, squares ):
        self.match = match
        self.ss_squares = squares
        self.ss_squares_data = []
    
    def compute( self ):
        self.matrix = []
        self.curve_row_mean = []
        self.curve_col_mean = []
        self.curve_local_rmsd = []
        
        # build profile
        for i in range( self.match.get_length() ):
            dp_util.Msg.out( "%sstep: %d of %d" %("\b" * 40, i, self.match.get_length()-1) )
            
            # superimposes 'cmp' on 'ref' minimizing RMSD of residue 'i'
            (ra1, ca1) = self.match.get_atoms( i )
            
            sup = Superimposer()
            sup.set_atoms( ra1, ca1 )
            sup.apply( self.match.cmp_model.get_atoms() )
            
            self.curve_local_rmsd.append( sup.rms )
            
            # Distance normalization
            if( NORMALIZE ):
                cm1 = self.center_of_mass(ra1) 
            
            # for each nucleotide
            row = []
            for j in range( self.match.get_length() ):
                dist_count = 0.0
                dist_sum = 0.0
                
                (ra2, ca2) = self.match.get_atoms( j )
                
                for k in X_( ra2 ):
                    dist = ra2[k] - ca2[k]
        
                    dist_count += 1.0
                    dist_sum += dist

                # append the average distance between all atoms of the given nucleotide
                # it gives how distant the ref and cmp nucleotides are when aligning the pivot nucleotide
                
                # Distance normalization
                if( NORMALIZE  and (i != j) ):
                    cm2 = self.center_of_mass(ra2)
                    
                    dist_norm = np.sqrt(np.sum((cm1-cm2) **2))
                    row.append( dist_sum / (dist_count * dist_norm) )
                else:
                    row.append( dist_sum / dist_count )
            
            self.matrix.append( row )
            self.curve_row_mean.append( np.mean( row ) )

        # full matrix
        self.matrix = np.array( self.matrix )
        
        # compute column mean
        for c in range( self.match.get_length() ):
            self.curve_col_mean.append( np.mean( self.matrix[:,c] ) )
            
        self.curve_row_mean = np.array( self.curve_row_mean )
        self.curve_col_mean = np.array( self.curve_col_mean )
        self.curve_local_rmsd = np.array( self.curve_local_rmsd )
        
        # get squares data
        self.compute_squares_data()
        dp_util.Msg.out( "%sdone\n" %("\b" * 40) )
    
    def center_of_mass(self, atoms):
        result = np.array( [0.0, 0.0, 0.0] )
        
        for atom in atoms:
            result += atom.get_vector().get_array()

        return( result / len(atoms) )

    def compute_squares_data(self):
        self.ss_squares_data = []
        
        # squares (r1, c1)-(r2,c2) => (a, b)-(c, d)
        for ( (a, b, c, d), txt, color) in self.ss_squares:
            width = c - a + 1
            height = d - b + 1

            size = width * height
            total = self.matrix[b:d+1,a:c+1].sum()
            avg = total / float(size)
            
            self.ss_squares_data.append( (a, b, c, d, width, height, txt, color, avg, total, size) )

    def matrix_save(self, fname):
        # write sequence text
        txt = "#DP 1.0\n"
        
        (rntxt, cntxt) = ("", "")
        (rrtxt, crtxt) = ("", "")
        for i in range( self.match.get_length() ):
            (rname, cname) = self.match.get_res_names( i )
            rntxt += rname
            cntxt += cname
            
            (rres, cres) = self.match.get_residues( i )
            rrtxt += "(%s:%s:'%s')" %(rres.get_parent().id, rres.get_id()[1], rres.get_resname())
            crtxt += "(%s:%s:'%s')" %(cres.get_parent().id, cres.get_id()[1], cres.get_resname())
                
        txt += "REF_PDB\t%s\nREF_MODEL\t%d\n" %(self.match.ref_pdb, self.match.ref_model_id)
        txt += "REF_MODEL_SEQUENCE\t%s\n" %rntxt
        txt += "REF_MODEL_RESIDUES\t%s\n" %rrtxt
        
        txt += "CMP_PDB\t%s\nCMP_MODEL\t%d\n" %(self.match.cmp_pdb, self.match.cmp_model_id)
        txt += "CMP_MODEL_SEQUENCE\t%s\n" %cntxt
        txt += "CMP_MODEL_RESIDUES\t%s\n" %crtxt
        
        # write profiles
        txt += "LOCAL_RMSD\t%s\n" %("\t".join( map( lambda x: "%.3f" %x, self.curve_local_rmsd ) ))
        txt += "ROW_MEANS\t%s\n" %("\t".join( map( lambda x: "%.3f" %x, self.curve_row_mean ) ))
        txt += "COL_MEANS\t%s\n" %("\t".join( map( lambda x: "%.3f" %x, self.curve_col_mean ) ))
        
        # write squares
        data_total = {}
        data_size = {}
        for (a, b, c, d, width, height, k, color, avg, total, size) in self.ss_squares_data:
            data_total[k] = data_total.get( k, 0.0 ) + total
            data_size[k] = data_size.get( k, 0 ) + size

        keys = data_total.keys()
        keys = list(keys)  # Convert dict_keys to a list
        keys.sort()
        for k in keys:
            txt += "SQUARE_VALUE\t'%s'\t%.3f\n" %( k, data_total[k]/float(data_size[k]) )
        
        # write matrix
        for i in range( self.match.get_length() ):
            txt += "ROW_%d\t%s\n" %(i, "\t".join( map( lambda x: "%.3f" %x, self.matrix[i] ) ))
            
        txt += "#eof"
        
        fo = open( fname, "w" )
        fo.write( txt )
        fo.close()
    
    def svg_save(self, fname):
        g = Graphics(self, fname, self.ss_squares_data )
        g.draw()

class Palette:
    def __init__(self, limit_down, steps_down, limit_up, steps_up):
        self.limit_down = limit_down
        self.steps_down = steps_down
        self.limit_up = limit_up
        self.steps_up = steps_up
        self.colors = []
        
        # From WHITE to YELLOW in 'self.steps_down' steps
        for i in range(self.steps_down, 0, -1):
            self.colors.append( svg.colorstr( 255, 255, int(float(i)/self.steps_down*255.0) ) )
    
        # From YELLOW to RED in 'self.steps_up' steps
        for i in range(self.steps_up, 0, -1):
            self.colors.append( svg.colorstr( 255, int(float(i)/self.steps_up*255.0), 0 ) )

    def get_colors(self, value):
        # choose the color
        if( value < self.limit_down ):
            ndx = min( int( value * float(self.steps_down) / self.limit_down ), self.steps_down - 1 )
        else:
            ndx = self.steps_down + min( int( (value - self.limit_down) * float(self.steps_up) / (self.limit_up-self.limit_down) ), self.steps_up-1 )
            
        return( self.colors[ndx] )

    
class Graphics:
    def __init__(self, dp, fname, squares):
        self.dp = dp
        self.fname = fname
        self.squares = squares
        
    def draw(self):
        self.prepare_stage()
        self.prepare_palette()
        
        self.draw_scale()
        self.draw_curves()
        self.draw_matrix()
        self.draw_secondary_structure()
        
        self.scene.write_svg( self.fname )

    def prepare_stage(self):
        self.scale_start = (SQR_SIDE, (17 + len(self.dp.matrix)) * SQR_SIDE)
        self.curve_start = (SQR_SIDE, (13 + len(self.dp.matrix)) * SQR_SIDE)
        self.matrix_start = (SQR_SIDE, (1 + len(self.dp.matrix)) * SQR_SIDE)
        
        self.scene = svg.Scene( (17 + len(self.dp.matrix)) * SQR_SIDE, (2 + len(self.dp.matrix)) * SQR_SIDE )
        
    def prepare_palette(self):
        m = copy.copy( self.dp.matrix )
        m = np.reshape( m, m.size )
        m.sort()
        
        # *** TODO ***: these values must be set in the config file
        #self.palette = Palette(m[m.size * LIMIT_DOWN], STEPS_DOWN, m.max(), STEPS_UP)
        if( NORMALIZE ):
            self.palette = Palette(1.5, STEPS_DOWN, 3, STEPS_UP)
        else:
            self.palette = Palette(30, STEPS_DOWN, 60, STEPS_UP)

    def draw_scale(self):
        self.coords_set( self.scale_start )
        
        x0 = 0
        x1 = self.palette.steps_down
        x2 = self.palette.steps_down + self.palette.steps_up
        
        # Color scale
        for i in range(x2):
            coords = self.coords( i * SQR_SIDE, 2 * SQR_SIDE )
            self.scene.add( svg.Rectangle( coords, SQR_SIDE, SQR_SIDE, self.palette.colors[i] ) )
        
        origin = self.coords( x0, 0 )
        self.scene.add( svg.Text( origin, "0.0" ) )

        origin = self.coords( x1 * SQR_SIDE, 0 )
        self.scene.add( svg.Text( origin, "%.2f" %self.palette.limit_down ) )
        
        origin = self.coords( (x2-1) * SQR_SIDE, 0 )
        self.scene.add( svg.Text( origin, "%.2f" %self.palette.limit_up ) )
        
    def draw_curves(self):
        self.coords_set( self.curve_start )
        
        x1 = len(self.dp.matrix)
        
        y_max_inc = max( self.dp.curve_row_mean.max(), self.dp.curve_col_mean.max(), self.dp.curve_local_rmsd.max() ) / 10.0
        y_max_inc = ((int(y_max_inc) / 5) + 1) * 5
        
        start = self.coords( 0, 0 )

        # Draw the XX axis
        end = self.coords( x1 * SQR_SIDE, 0 )
        self.scene.add( svg.Line( start, end ) )
        
        # draw both YY axis
        end = self.coords( 0, 10 * SQR_SIDE )
        self.scene.add( svg.Line( start, end ) )
        
        # draw scales in 5 steps
        for i in range( 0, 6 ):
            origin = self.coords( -SQR_SIDE, i * 2 * SQR_SIDE )
            self.scene.add( svg.Text( origin, "%2d" %round(y_max_inc * 2 * float(i)) ), 10 )

        # draw horizontal scale
        for r in X_(self.dp.matrix):
            (rname, cname) = self.dp.match.get_res_names( r )
            origin = self.coords( r * SQR_SIDE, -SQR_SIDE )
            self.scene.add( svg.Text( origin, cname ) )

        # draw curves
        self.draw_curves_aux( self.dp.curve_row_mean, y_max_inc , x1, COLOR_ROW_MEAN )
        self.draw_curves_aux( self.dp.curve_col_mean, y_max_inc , x1, COLOR_COL_MEAN )
        self.draw_curves_aux( self.dp.curve_local_rmsd, y_max_inc , x1, COLOR_LOCAL_RMSD )
        
    def draw_curves_aux(self, points, y_max_inc, x1, stroke_color):
        point_0 = self.coords( SQR_SIDE / 2, (points[0]/y_max_inc) * SQR_SIDE )
        
        for i in range( 1, x1 ):
            point_1 = self.coords( i * SQR_SIDE + SQR_SIDE / 2, (points[i]/y_max_inc) * SQR_SIDE )

            self.scene.add( svg.Line( point_0, point_1, stroke_color=stroke_color, stroke_width=1 ) )
            
            point_0 = point_1

    def draw_matrix(self):
        self.coords_set( self.matrix_start )
                
        # draw matrix and write sequence
        for r in X_(self.dp.matrix):
            (rname, cname) = self.dp.match.get_res_names( r )
            
            origin = self.coords( -SQR_SIDE, r * SQR_SIDE )
            self.scene.add( svg.Text( origin, rname ) )
            #self.scene.add( svg.Text( origin, r ) )
            
            origin = self.coords( r * SQR_SIDE, -SQR_SIDE )
            self.scene.add( svg.Text( origin, cname ) )
            #self.scene.add( svg.Text( origin, r ) )
            
            for c in X_(self.dp.matrix):
                origin = self.coords( c * SQR_SIDE, (r+1) * SQR_SIDE )
                color = self.palette.get_colors( self.dp.matrix[r,c] )
                self.scene.add( svg.Rectangle( origin, SQR_SIDE, SQR_SIDE, color, color, 0) )
       
    def draw_secondary_structure(self):
        # squares (r1, c1)-(r2,c2) => (a, b)-(c, d)
        for ( a, b, c, d, width, height, txt, color, avg, dummy, dummy ) in self.squares:
            origin = self.coords( a*SQR_SIDE+1, (d+1)*SQR_SIDE-1 )
            self.scene.add( svg.Rectangle( origin, (height * SQR_SIDE)-2, (width * SQR_SIDE)-2, stroke_color=color, stroke_width=2), 5 )
            
            origin = self.coords( a*SQR_SIDE+2, d*SQR_SIDE )
            self.scene.add( svg.Text( origin, str(txt), size=8 ), 5 )
            
            origin = self.coords( a*SQR_SIDE+2, b*SQR_SIDE+3 )
            self.scene.add( svg.Text( origin, "%.2f" %avg, size=8 ), 5 )
                    
    def coords_set(self, start):
        self.coords_current_start = start
        
    def coords(self, x, y):
        return( (self.coords_current_start[0] + x, self.coords_current_start[1] - y) )