# --------------------------------------------------------------------
# dp_2d.py
#
# Secondary structure model class.
#
# history:
#    20090825 - 1.0.0 - JAC - first version
#    20090825 - 1.0.1 - JAC - minor bug corrections
# --------------------------------------------------------------------

import dp_util

class SecondaryStructure:
    def __init__(self, helices=[], loops=[]):
        self.helices = {}
        for h in helices:
            self.helices[h.name] = h

        self.loops = {}
        for l in loops:
            self.loops[l.name] = l

    def square_loop(self, name, square_name=""):
        square_name = (square_name == "" and name) or square_name

        result = []
        l = self.loops[name]
        
        (a, b, c, d) = (l.i, l.i, l.i + l.n - 1, l.i + l.n - 1)
        
        result.append( (SecondaryStructure.get_upper(a, b, c, d), square_name, "#ff00ff") )
        
        for sq in result:
            SecondaryStructure.check_coords( sq[0], "loop %s" %name, True )
        
        return( result )
              
    def square_helix(self, name, square_name=None, upper=None):
        square_name = (square_name == "" and name) or square_name
        
        result = []
        h = self.helices[name]
        
        (a, b, c, d) = (h.i, h.j, h.i + h.ni - 1, h.j + h.nj - 1)
        
        if( upper==None or upper ):
            result.append( (SecondaryStructure.get_upper(a, b, c, d), square_name, "#ff00ff") )
        
        if( upper==None or not upper ):
            result.append( (SecondaryStructure.get_lower(a, b, c, d), square_name, "#ff00ff") )
        
        for sq in result:
            SecondaryStructure.check_coords( sq[0], "helix %s" %name )
        
        return( result )

    def square_hh(self, name1, name2, square_name=None, upper=None):
        square_name = (square_name == "" and "%s x %s" %(name1, name2)) or square_name
        
        result = []
        temp = []
        
        h1 = self.helices[name1]
        h2 = self.helices[name2]
        
        temp.append( (h1.i, h2.i, h1.i + h1.ni - 1, h2.i + h2.ni - 1) )    # strand h1a, h2a
        temp.append( (h1.i, h2.j, h1.i + h1.ni - 1, h2.j + h2.nj - 1) )    # strand h1a, h2b
        temp.append( (h1.j, h2.i, h1.j + h1.nj - 1, h2.i + h2.ni - 1) )    # strand h1b, h2a
        temp.append( (h1.j, h2.j, h1.j + h1.nj - 1, h2.j + h2.nj - 1) )    # strand h1b, h2b

        for t in temp:
            if( upper==None or upper ):
                result.append( (SecondaryStructure.get_upper(t[0], t[1], t[2], t[3]), square_name, "#0000ff") )
            
            if( upper==None or not upper ):
                result.append( (SecondaryStructure.get_lower(t[0], t[1], t[2], t[3]), square_name, "#0000ff") )

        for sq in result:
            SecondaryStructure.check_coords( sq[0], "comparing helices %s and %s" %(name1, name2) )

        return( result )

    def square_ll(self, name1, name2, square_name=None, upper=None):
        square_name = (square_name == "" and "%s x %s" %(name1, name2)) or square_name
                
        result = []
        
        l1 = self.loops[name1]
        l2 = self.loops[name2]
        
        (a, b, c, d) = (l1.i, l2.i, l1.i + l1.n - 1, l2.i + l2.n - 1)
        
        if( upper==None or upper ):
            result.append( (SecondaryStructure.get_upper(a, b, c, d), square_name, "#0000ff") )
        
        if( upper==None or not upper ):
            result.append( (SecondaryStructure.get_lower(a, b, c, d), square_name, "#0000ff") )
        
        for sq in result:
            SecondaryStructure.check_coords( sq[0], "comparing loops %s and %s" %(name1, name2) )
        
        return( result )

    def square_hl(self, name1, name2, square_name=None, upper=None):
        square_name = (square_name == "" and "%s x %s" %(name1, name2)) or square_name
        
        result = []
        temp = []
        
        h = self.helices[name1]
        l = self.loops[name2]
        
        temp.append( (l.i, h.i, l.i + l.n - 1, h.i + h.ni - 1) )    # strand h1a, h2a
        temp.append( (l.i, h.j, l.i + l.n - 1, h.j + h.nj - 1) )    # strand h1a, h2b
        
        for t in temp:
            if( upper==None or upper ):
                result.append( (SecondaryStructure.get_upper(t[0], t[1], t[2], t[3]), square_name, "#0000ff") )
            
            if( upper==None or not upper ):
                result.append( (SecondaryStructure.get_lower(t[0], t[1], t[2], t[3]), square_name, "#0000ff") )

        for sq in result:
            SecondaryStructure.check_coords( sq[0], "comparing helix %s and loop %s" %(name1, name2) )

        return( result )

    def square_lh(self, name1, name2, square_name=None, upper=None):
        return( self.square( name2, name1, square_name, upper ) )
    
    def get_upper(a, b, c, d):
        if( a > b ):
            return( b, a, d, c )
        else:
            return( a, b, c, d )

    def get_lower(a, b, c, d):
        if( a < b ):
            return( b, a, d, c )
        else:
            return( a, b, c, d )
        
    def check_coords(sq, err_str, single_strand=False):
        (a, b, c, d) = sq
        
        if( not single_strand ):
            # its XOR since the conditions must be both true or false
            if( bool(a>b) ^ bool(c>d) ^ bool(a>d) ^ bool(c>b) ):
                dp_util.Msg.fatal( "Bad cross over in %s\n(%d, %d) - (%d, %d)" %(err_str, a, b, c, d) )

        if( (a>c) or (b>d) ):
            dp_util.Msg.fatal( "Inverted Square in %s\n(%d, %d) - (%d, %d)" %(err_str, a, b, c, d) )

    
    check_coords = staticmethod(check_coords)
    get_upper = staticmethod(get_upper)
    get_lower = staticmethod(get_lower)
        
        
class Helix:
    def __init__(self, name, i, ni, j, nj):
        self.name = name
        self.i = i
        self.ni = ni
        self.j = j
        self.nj = nj
        
class Loop:
    def __init__(self, name, i, n):
        self.name = name
        self.i = i
        self.n = n