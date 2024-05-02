# --------------------------------------------------------------------
# svg.py
#
# SVG creator classes.
#
# To create an SVG class start with a 'Scene' object and new items like
# Circles, Lines, Rectangles and Text. When done just call 'write_svg'  
#
# history:
#    20090825 - 1.0.0 - JAC - first version
# --------------------------------------------------------------------
                
class Scene:
    def __init__(self, height=400, width=400 ):
        self.items = {}
        self.height = height
        self.width = width
        return

    def add(self,item,key=0):
        if key not in self.items:
            self.items[key] = []
            
        self.items[key].append(item)

    def xml(self):
        tag_svg = Tag( "svg" )
        tag_svg.add_attribute( "xmlns", "http://www.w3.org/2000/svg" )
        tag_svg.add_attribute( "xmlns:xlink", "http://www.w3.org/1999/xlink" )
        tag_svg.add_attribute( "height", self.height )
        tag_svg.add_attribute( "width", self.width )

        tag_g = Tag( "g" )
        tag_g.add_attribute( "style", "fill-opacity:1.0; stroke:black; stroke-width:1;" )
        
        tag_svg.add_child( tag_g )
        
        keys = self.items.keys()
        keys = list(keys)  # Convert dict_keys to a list
        keys.sort()
        
        for key in keys:
            for item in self.items[key]:
                tag_g.add_child( item.get_tags() )

        var = ""
        var += tag_svg.xml()
        
        return var

    def write_svg(self, filename=None, height=None, width=None ):
        if( height != None ):
            self.height = height
            
        if( width != None ):
            self.width = width
        
        if( not filename.endswith( ".svg" ) ):
            filename += ".svg"
            
        file = open( filename, 'w' )
        file.write( self.xml() )
        file.close()

class Line:
    def __init__(self, start, end, stroke_color="#000000", stroke_width=1):
        self.start = start
        self.end = end
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
    
    def get_tags(self):
        tag = Tag( "line" )
        tag.add_attribute( "x1", self.start[0] )
        tag.add_attribute( "y1", self.start[1] )
        tag.add_attribute( "x2", self.end[0] )
        tag.add_attribute( "y2", self.end[1] )
        tag.add_attribute( "stroke", self.stroke_color )
        tag.add_attribute( "stroke-width", self.stroke_width )
        
        return( tag )        
        

class Circle:
    def __init__(self, center, radius, color):
        self.center = center #xy tuple
        self.radius = radius #xy tuple
        self.color = color   #html color

    def get_tags(self):
        tag = Tag( "circle" )
        tag.add_attribute( "cx", self.center[0] )
        tag.add_attribute( "cy", self.center[1] )
        tag.add_attribute( "r", self.radius )
        tag.add_attribute( "style", "fill:%s;" %self.color )
        
        return( tag )        

class Rectangle:
    def __init__(self,origin,height,width,fill="none",stroke_color="#000000",stroke_width=1):
        self.origin = origin
        self.height = height
        self.width = width
        self.fill = fill
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
    
    def get_tags(self):
        tag = Tag( "rect" )
        tag.add_attribute( "x", self.origin[0] )
        tag.add_attribute( "y", self.origin[1] )
        tag.add_attribute( "height", self.height )
        tag.add_attribute( "width", self.width )
        tag.add_attribute( "fill", self.fill )
        tag.add_attribute( "stroke", self.stroke_color )
        tag.add_attribute( "stroke-width", self.stroke_width )
        
        return( tag )

class Text:
    def __init__(self,origin,text,size=10,color="#000000"):
        self.origin = origin
        self.text = text
        self.size = size
        self.color = color

    def get_tags(self):
        tag = Tag( "text" )
        tag.add_attribute( "x", self.origin[0] )
        tag.add_attribute( "y", self.origin[1] )
        tag.add_attribute( "font-size", self.size )
        tag.add_attribute( "font-weight", "normal" )
        tag.add_attribute( "font-family", "Verdana" )
        tag.add_attribute( "stroke", self.color )
        tag.add_cdata( self.text )
        
        return( tag )      

class Tag:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.attributes = []
        self.cdata = ""
        
    def add_child(self, tag):
        self.children.append( tag )
    
    def add_attribute(self, name, value):
        self.attributes.append( (name, value) )
        
    def add_cdata(self, cdata):
        self.cdata = cdata
        
    def xml(self):
        child_str = self.cdata
        attr_str = ""
        
        for child in self.children:
            child_str += child.xml()
        
        for attr in self.attributes:
            if( attr != None ):
                attr_str += " %s=\"%s\"" %(attr[0], str(attr[1]))
        
        if( child_str == "" ):
            result = "<%s %s/>\n" %(self.name, attr_str)
        else:
            result = "<%s %s>\n%s</%s>\n" %(self.name, attr_str, child_str, self.name)
        
        return( result )
    
def colorstr( r, g, b ):
    return "#%02x%02x%02x" % (r,g,b)

def test():
    scene = Scene()
    scene.add(Rectangle((100,100),200,200,"#00ffff"))
    scene.add(Line((200,200),(200,300)))
    scene.add(Line((200,200),(300,200)))
    scene.add(Line((200,200),(100,200)))
    scene.add(Line((200,200),(200,100)))
    scene.add(Circle((0,0),30,colorstr(0,0,255)))
    scene.add(Text((50,50),"Testing SVG"))
    scene.write_svg('i:\\testx.svg')

if __name__ == '__main__':
    test()
