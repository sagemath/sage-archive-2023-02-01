"""nodoctest
Package for drawing graphs related to elliptic curves, etc.

EXAMPLES:

g=Graph('ex2', (-3,-5),(3,5))
g.color(0.7)
g.grid()
g.color(0.4)
g.axes()
g.color(0,0,0.7)
g.elliptic_curve([0,0,1,-1,0])
g.dvi()

"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# standard Python modules
import math, os, pickle

# Tools used for viewing files.  This list should probably be in a
# central location, not in this file.

tools = {"PS_VIEWER":"gv",\
         "LATEX":"latex",\
         "DVIPS":"dvips", \
         "GRAPHS_DIR":"graphs", \
         "CONVERT":"convert", \
         "PS2PDF": "ps2pdf"}

class Point:
    """

    """
    def __init__(self, x, y=None):
        if isinstance(x, tuple):
            assert len(x) == 2
            x, y = x
        elif isinstance(x, Point):
            x, y = x.tuple()
        self.__x = float(x)
        self.__y = float(y)

    def __getitem__(self, n):
        if n == 0:
            return self.__x
        elif n == 1:
            return self.__y
        else:
            raise IndexError

    def __add__(self, other):
        return Point(self.__x + other.__x, self.__y + other.__y)

    def __cmp__(self, other):
        if not isinstance(other, Point):
            return -1
        if self.__x < other.__x:
            return -1
        elif self.__x > other.__x:
            return 1
        elif self.__y < other.__y:
            return -1
        elif self.__y > other.__y:
            return +1
        else:
            return 0

    def __repr__(self):
        return "(%.3f,%.3f)"%(self.__x,self.__y)

    def __rmul__(self, scalar):
        return Point(scalar*self.__x, scalar*self.__y)

    def __sub__(self, other):
        return Point(self.__x - other.__x, self.__y - other.__y)

    def d(self, other):
        return math.sqrt((self.__x-other.__x)**2 + (self.__y - other.__y)**2)

    def tuple(self):
        return (self.__x, self.__y)

    def x(self):
        return self.__x

    def y(self):
        return self.__y





class Command:

    def __init__(self, x):
        self.__repr = x

    def __repr__(self):
        return self.rep

    def expand(self):
        return self.exp


class Color(Command):

    def __init__(self, red, green=None, blue=None, colorname="mycolor"):
        if green == None and blue == None:
            green, blue = red, red
        if isinstance(red, Color):
            if red.__rgb == None:
                self.rep = red.rep
                self.exp = red.exp
                return
            red, green, blue = red.__rgb
        elif isinstance(red, tuple):
            red, green, blue = red
        elif isinstance(red, str):
            self.rep = "Color "+red
            self.exp = "\\%s"%red
            self.__rgb = None
            return
        self.__rgb = (float(red), float(green), float(blue))
        self.rep = "RGB color " + str(self.__rgb)
        self.exp = "\\newrgbcolor{%s}{%.2f %.2f %.2f}"%(
            colorname, self.__rgb[0], self.__rgb[1], self.__rgb[2])
        if colorname=="mycolor":
            self.exp += "\\%s\n"%colorname

    def rgb(self):
        return self.__rgb

class Line(Command):

    def __init__(self, start, stop, style="-"):
        self.rep = "Line from %s to %s"%(start, stop)
        self.exp = "\\psline[linecolor=mycolor]{%s}%s%s"%(style,start,stop)


class Grid(Command):

    def __init__(self, subgriddiv=1, gridlabelcolor=Color(0)):
        self.rep = "Grid divided into %s subdivisions"%\
                   (subgriddiv)
        self.exp = Color(gridlabelcolor,colorname="glc").expand() + \
                   "\\psgrid[gridcolor=mycolor, subgriddiv=%s, gridlabelcolor=glc]"%(subgriddiv)

class Set(Command):

    def __init__(self, property, value):
        self.rep = "%s = %s"%(property, value)
        self.exp = "\\psset{%s=%s}"%(property, value)


class Text(Command):

    def __init__(self, position, value, refh="l", refv="b"):
        """
        refh = l (left), r (right), c (center)
        refv = t (top), b (bottom), B (baseline)
        """
        if refh=="c":
            refh=""
        self.rep = '"%s" at position %s'%(value, position)
        self.exp = '\\rput[%s%s]%s{%s}'%(refh, refv, position, value)

class Curve(Command):

    def __init__(self, points, max_points=0):
        self.rep = "Curve defined by %s points"%len(points)
        if len(points) == 0:
            self.exp = ""
        else:
            self.exp = "\\pscurve[linecolor=mycolor]"
            if max_points > 0:
                n = len(points)/float(max_points)
                if n > 1:
                    n = int(n+1)
                    points = [points[n*k] for k in range(len(points)/n)]
            for x in points:
                self.exp += str(x)

class Function(Command):
    """
    Draw
    """
    def __init__(self, f, xmin, xmax, num_points=None, smooth=True,
                 shade_below=False, shade_color=0, line_color=0,
                 y_cutoff=None, dot=0, y_scale=1.0, progress=False):
        self.rep = "Graph of a function."
        if num_points == None:
            num_points = 100
        s = float(xmax - xmin)
        interval = s / num_points
        if progress:
            self.points = []
            for i in xrange(num_points + 1):
                print i
                self.points.append((xmin + i*interval, y_scale*f(xmin + i * interval)))
        else:
            self.points = [(xmin + i*interval, y_scale*f(xmin + i * interval)) for
                           i in xrange(num_points+1)]
        if y_cutoff != None:
            self.points = [(x,y) for x,y in self.points if abs(y) <= y_cutoff]
        points = self.points
        self.exp = ""
        if shade_below:
            closed = [(xmin,0)] + points + [(xmax,0)] + [(xmin,0)]
            self.closed = closed
            self.exp += Color(shade_color).exp + "\n"
            self.exp += Polygon(closed, filled=True).exp + "\n"
        self.exp += Color(line_color).exp + "\n"
        if smooth:
            self.exp += Curve(points).exp + "\n"
        else:
            self.exp += Path(points).exp + "\n"
        if dot>0:
            for P in points:
                self.exp += Disc(P,dot).exp + "\n"


class Path(Command):

    def __init__(self, points, style="-", max_points=0):
        self.rep = "Path defined by %s points"%len(points)
        if len(points) == 0:
            self.exp = ""
        else:
            self.exp = ""
            if max_points > 0:
                n = len(points)/float(max_points)
                if n > 1:
                    n = int(n+1)
                    points = [points[n*k] for k in range(len(points)/n)]

            for i in range(len(points)-1):
                P = points[i]
                Q = points[i+1]
                self.exp += "\n\\psline[linecolor=mycolor]{%s}(%s,%s)(%s,%s)"%(
                    style, P[0], P[1], Q[0], Q[1])

class Disc(Command):

    def __init__(self, center, radius):
        if not isinstance(center, Point):
            center = Point(center)
        if not isinstance(radius,float):
            radius = float(radius)
        self.rep = "Solid dot at %s of radius %s"%(center, radius)
        self.exp = "\\pscircle*[linecolor=mycolor]%s{%.3f}"%(center, radius)

Dot = Disc   # for backward compatibility with some code I wrote.

class Circle(Command):

    def __init__(self, center, radius):
        self.rep = "Circle at %s of radius %s"%(center, radius)
        self.exp = "\\pscircle[linecolor=mycolor]%s{%s}"%(center, radius)


class Polygon(Command):

    def __init__(self, points, filled=True):
        """
        INPUT:
            points -- a list of points
            filled -- True or False (default is True)
        """
        self.pts = [Point(x) for x in points]
        self.rep = "Polygon with vertices"%self.pts
        X = "".join([str(x) for x in self.pts])
        if filled:
            fill = "*"
        else:
            fill = ""
        self.exp = "\\pspolygon%s[linecolor=mycolor]"%fill + X


class Axes(Command):

    def __init__(self, graph, x="$x$", y="$y$", offset=0.2):
        xmin, ymin = graph.lower_left().tuple()
        xmax, ymax = graph.upper_right().tuple()
        unit = graph.unit()
        self.rep = "%s-%s axes"%(x,y)
        self.exp = "%s\n%s\n%s\n%s\n%s\n%s\n"%\
                   (Line(Point(xmin,0),Point(xmax,0),"->").exp, \
                        Line(Point(0,ymin), Point(0,ymax), "->").exp, \
                        Color(0).exp,
                        Text(Point(xmax+offset,0),"{\\Large %s}"%x).exp, \
                        Text(Point(0,ymax+offset),"{\\Large %s}"%y).exp, \
                        Color(graph.current_color).exp)

class Ticks(Command):
    """
    Put a vertical tick mark and label at each position
    on the x-axis between xmin and xmax that when scaled
    is an integer multiple of xmult.  Note that
    xmult must be an integer.  Same with y axis.
    """
    def __init__(self, \
                 xmin, xmax, xscale, xmult, xlen, above, \
                 ymin, ymax, yscale, ymult, ylen, right):
        self.rep="Axes tick marks"
        self.exp = ""
        assert isinstance(xmult,int)

        n = int(xscale*(xmax - xmin)/float(xmult))+1  # number of tick marks
        xstart = int(math.ceil(xmin*xscale))
        xstart += xmult*(xstart/xmult) - xstart
        assert xstart%xmult == 0
        if above:
            ypos = xlen
        else:
            ypos = -3*xlen
        for i in range(n):
            scaled = xstart + i*xmult
            if scaled == 0:
                continue
            unscaled = scaled/float(xscale)
            self.exp += "%s\n"%(Line((unscaled,-xlen/2.0),(unscaled,xlen/2.0)).exp)
            self.exp += "%s\n"%(Text((unscaled,ypos),"%s"%scaled, refh="c").exp)

        # Same but for y axis
        n = yscale*(ymax - ymin)/float(ymult) + 1 # number of tick marks
        ystart = int(math.ceil(ymin*yscale))
        ystart += ymult*(ystart/ymult) - ystart
        assert ystart%ymult == 0
        if right:
            ypos = 3*ylen
        else:
            ypos = -3*ylen
        for i in range(int(n)):
            scaled = ystart + i*ymult
            if scaled == 0:
                continue
            unscaled = scaled/float(yscale)
            self.exp += "%s\n"%(Line((-ylen/2.0, unscaled),(ylen/2.0, unscaled)).exp)
            self.exp += "%s\n"%(Text((ypos, unscaled-ylen),"%s"%scaled, refh="c").exp)



class EllipticCurve(Command):

    def __init__(self, graph, a1,a2,a3,a4,a6, res=15):
        assert isinstance(graph, Graph)
        assert isinstance(a1,(int,float))
        assert isinstance(a2,(int,float))
        assert isinstance(a3,(int,float))
        assert isinstance(a4,(int,float))
        assert isinstance(a6,(int,float))
        assert isinstance(res,(int,float))

        self.rep = "Elliptic curve with invariants [%s,%s,%s,%s,%s]"%\
                   (a1,a2,a3,a4,a6)
        self.exp = ""
        unit = graph.unit()
        delta = unit/float(res)
        xmin, ymin = graph.lower_left().tuple()
        xmax, ymax = graph.upper_right().tuple()
        def y(x):
            try:
                r = math.sqrt((a1*x+a3)**2 + 4.0*(x**3+a2*x**2+a4*x+a6))
            except ValueError:
                return (None,None)
            return ( Point(x,(-(a1*x+a3)-r)/2.0), Point(x,(-(a1*x+a3)+r)/2.0) )

        # 1. Compute the non-compact component
        xx = xmax
        x0 = xx
        P,_ = y(xx)
        P0 = P
        C1 = []
        while P != None:
            C1.append(P)
            xx -= delta
            Q, R = y(xx)
            if Q == None:
                x0 = xx
                break
            P0 = P
            if P.d(Q) < P.d(R):
                P = Q
            else:
                P = R
        P = P0
        while xx <= xmax:
            xx += delta
            Q, R = y(xx)
            if Q == None:
                continue
            if P.d(Q) < P.d(R) and not (Q in C1):
                P = Q
            else:
                P = R
            C1.append(P)
        #end while
        self.exp += Curve(C1).exp
        # Now draw the second compact component, if it exists:
        xx = x0
        Q, R = y(xx)
        while xx >= xmin and Q == None:
            xx -= delta
            Q, R = y(xx)
        if xx == xmin or Q==None:  # no compact component
            return
        C2 = []
        P = Q
        first = P
        while P != None:
            C2.append(P)
            xx -= delta
            Q, R = y(xx)
            if Q == None:
                break
            P0 = P
            if P.d(Q) < P.d(R):
                P = Q
            else:
                P = R
            C2.append(P)
        #end while
        P = P0
        while P != None:
            xx += delta
            Q, R = y(xx)
            if Q == None:
                break
            if P.d(Q) < P.d(R) and not (Q in C2):
                P = Q
            else:
                P = R
            C2.append(P)
        #end while
        C2.append(first)
        self.exp += Curve(C2).exp



def load_graph(filename):

    name = "%s/%s.pickle"%(tools["GRAPHS_DIR"], filename)
    if os.path.exists(name):
        return pickle.load(open(name,"r"))
    raise RuntimeError, "File %s not found."%name

class Graph:
    def __init__(self, name, lower_left=Point(0,0),
                 upper_right=Point(10,10), unit=None):
        self.__name = name
        if not isinstance(lower_left, Point):
            lower_left = Point(lower_left)
        if not isinstance(upper_right, Point):
            upper_right = Point(upper_right)
        if unit==None:
            span = upper_right - lower_left
            unit = 25.0/max(span.x(), span.y())
        if not os.path.exists(tools["GRAPHS_DIR"]):
            os.mkdir(tools["GRAPHS_DIR"])
        self.__reset()
        self.current_color = Color(0)
        self.__pen = None
        self.__grid = None
        self.__unit = unit
        self.__lower_left = Point(lower_left)
        self.__upper_right = Point(upper_right)
        self.__commands = []

    def __delitem__(self, n):
        del self.__commands[n]

    def __getitem__(self, n):
        return self.__commands[n]

    def __repr__(self):
        s = "Graph '%s':\n"%self.__name
        i = 1
        for x in self.commands():
            s += str(i) + ":\t" + str(x) + "\n"
        return s

    def __reset(self):
        self.__done = []

    def axes(self, x="$x$", y="$y$", offset=0.2):
        self.append(Axes(self, x, y, offset=offset))

    def append(self, cmd):
        self.__reset()
        self.__commands.append(cmd)

    def circle(self, center, radius):
        self.append(Circle(center, radius))

    def color(self, color, blue=None, green=None):
        if isinstance(color, Color):
            self.append(color)
        else:
            if blue == None or green == None:
                blue = color
                green = color
            self.append(Color(color, blue, green))
        self.current_color = color

    def commands(self):
        return self.__commands

    def copy(self, name):
        x = Graph(name, self.__lower_left, self.__upper_right)
        for k in self.__dict__.keys():
            x.__dict__[k] = self.__dict__[k]
        x.__name = name
        x.__commands = list(self.__commands)
        return x

    def curve(self, points, max_points=0):
        self.append(Curve(points, max_points=max_points))

    def dot(self, center, radius=None):
        if radius==None:
            radius = self.__unit/10.0
        self.append(Dot(center, radius))

    def disc(self, center, radius=None):
        if radius==None:
            radius = self.__unit/10.0
        self.append(Disc(center, radius))

    def dvi(self):
        if "dvi" in self.__done:
            return
        self.tex()
        os.system('cd %s; %s "%s"'%\
                  (tools["GRAPHS_DIR"], tools["LATEX"], self.name()))
        assert os.path.exists("graphs/"+self.name()+".dvi"),\
               "Failed to create dvi file."
        self.__done.append("dvi")

    def elliptic_curve(self, a1=0,a2=0,a3=0,a4=0,a6=0, res=15.0):
        if isinstance(a1,list):
            a2=a1[1]; a3=a1[2]; a4=a1[3]; a6=a1[4]; a1=a1[0]
        self.append(EllipticCurve(self, a1,a2,a3,a4,a6, res))

    def emacs(self):
        self.tex()
        os.system('cd %s; e "%s.tex"&'%(tools["GRAPHS_DIR"], self.name()))

    def eps(self):
        if "eps" in self.__done:
            return
        self.ps()
        os.system('cd %s; %s -crop 0x0  "%s.ps" "%s.eps"'%\
                  (tools["GRAPHS_DIR"], tools["CONVERT"], self.name(),self.name()))
        self.__done.append("eps")

    def linestyle(self, style="solid"):
        """
        """
        self.append(Set("linestyle", style))

    def fillstyle(self, style="solid"):
        """
        Valid styles are
            none, solid, vlines, vlines*, hlines, hlines*,
            crosshatch and crosshatch*
        The default is solid.

        vlines, hlines and crosshatch draw a pattern of lines,
        according to the four parameters list below that are prefixed
        with hatch.  The * versions also fill the background, as
        in the solid style.

        The following are also methods that are relevant for
        hatch fills:

        hatchwidth = dim (default .8pt)
           Width of lines
        hatchsep=dim (default 4pt)
           Width of space between lines
        hatchcolor=color (default black)
           Color of lines
        hatchangle=rot (default: 45)
           Rotation of the lines, in degrees.  For example, if hatchangle is
           set to 45, the vlines style draws lines that run northwest-southeast,
           and the hlines style draws lines that run southwest-northeast,
           and the crosshatch style draws both together.
        """
        self.append(Set("fillstyle", style))

    def fillcolor(self, color):
        self.append(Color(color))
        self.append(Set("fillcolor", "mycolor"))

    def hatchwidth(self, dim="0.8pt"):
        self.append(Set("hatchwidth",dim))

    def hatchsep(self, dim="4pt"):
        self.append(Set("hatchsep",dim))

    def hatchcolor(self, color):
        self.append(Color(color))
        self.append(Set("hatchcolor","mycolor"))

    def hatchangle(self, rot=45):
        self.append(Set("hatchangle",rot))

    def gif(self, density=150):
        if "gif" in self.__done:
            return
        self.ps()
        print "Converting to gif"
        os.system('cd %s; %s -crop 0x0 -density %sx%s "%s.ps" "%s.gif"'%\
                  (tools["GRAPHS_DIR"], tools["CONVERT"], density, \
                   density, self.name(),self.name()))
        print "Done converting"
        self.__done.append("gif")

    def grid(self, subgriddiv=1, gridlabelcolor=Color(0)):
        self.append(Grid(subgriddiv, gridlabelcolor=gridlabelcolor))

    def jpg(self, density=150, quality=100):
        if "jpg" in self.__done:
            return
        self.ps()
        print "Converting to jpg"
        os.system('cd %s; %s -quality %s -crop 0x0 -density %sx%s "%s.ps" "%s.jpg"'%\
                  (tools["GRAPHS_DIR"], tools["CONVERT"], quality,
                   density, density, self.name(),self.name()))
        print "Done converting"
        self.__done.append("jpg")

    def line(self, start_point, stop_point, style="-"):
        self.append(Line(start_point, stop_point, style))

    def line_width(self, width):
        self.append(Set("linewidth",width))

    def function(self, f, xmin=None, xmax=None, num_points=None, smooth=True,
                 include_in_graph = True, shade_below=False, shade_color=0,
                 y_cutoff=None, dot=0, y_scale=1.0, progress=False):
        """
        ...
        Points with y coordinate bigger in abs value than y_cutoff are ignored
        (unless y_cutoff is not specified).
        """
        if xmin == None:
            xmin = self.__lower_left[0]
        if xmax == None:
            xmax = self.__upper_right[0]
        F = Function(f, xmin, xmax, num_points, smooth,
                     shade_below, shade_color, line_color=self.current_color,
                     y_cutoff = y_cutoff, dot=dot, y_scale=y_scale, progress=progress)
        if include_in_graph:
            self.append(F)
        return F

    def lower_left(self):
        return self.__lower_left

    def name(self):
        return self.__name

    def path(self, points, style="-", max_points=0):
        self.append(Path(points, style, max_points))

    def polygon(self, points, filled=True):
        """
        INPUT:
            points -- a list of points
            filled -- True or False (default is True)
        """
        self.append(Polygon(points, filled=filled))

    def png(self, density=150):
        if "png" in self.__done:
            return
        self.ps()
        print "Converting to png"
        os.system('cd %s; %s -crop 0x0 -density %sx%s "%s.ps" "%s.png"'%\
                  (tools["GRAPHS_DIR"], tools["CONVERT"], density,
                   density, self.name(),self.name()))
        print "Done converting"
        self.__done.append("png")

    def pdf(self):
        if "pdf" in self.__done:
            return
        self.ps()
        os.system('cd %s; %s "%s.ps" "%s.pdf"'%\
                  (tools["GRAPHS_DIR"], tools["PS2PDF"], self.name(),self.name()))
        self.__done.append("pdf")

    def ps(self):
        if "ps" in self.__done:
            return
        self.dvi()
        os.system('cd %s; %s -Ppdf -f < "%s.dvi" > "%s.ps"'%\
                  (tools["GRAPHS_DIR"], tools["DVIPS"], self.name(), self.name()))
        self.__done.append("ps")

    def save(self):
        fname = "graphs/%s.pickle"%self.__name
        pickle.dump(self, open(fname,"w"), 1)

    def tex(self):
        if "tex" in self.__done:
            return
        f = open("graphs/"+self.name()+".tex","w")
        f.write("\\documentclass{article}\n")
        f.write("\\pagestyle{empty}\n")
        f.write("\\thispagestyle{empty}\n")
        f.write("\\usepackage{amsmath}\n")
        f.write("\\usepackage{fullpage}\n")
        f.write("\\usepackage{pstricks}\n")
        f.write("\\begin{document}\n")
        f.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" + \
                "%%%%%%%%%%%%%%%%%%%%%%\n")
        f.write("%% Graph: %s\n"%self.name())
        f.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" + \
                "%%%%%%%%%%%%%%%%%%%%%%\n")
        f.write("\\psset{unit=%s}\n"%self.__unit)
        pspic = "%s_pic"%self.name()
        f.write("\\include{%s}\n"%pspic)
        g = open("graphs/" + pspic + ".tex", "w")
        g.write("%%\\psset{unit=%s}\n"%self.__unit)
        g.write("\\pspicture%s%s\n"%(self.__lower_left, self.__upper_right))
        g.write(Color(0,0,0).expand())
        for c in self.commands():
            g.write("%" + "%s\n%s\n"%(c,c.expand()))
        g.write("\\endpspicture\n")
        f.write("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" + \
                "%%%%%%%%%%%%%%%%%%%%%%\n")
        f.write("\\end{document}\n")
        self.__done.append("tex")

    def text(self, position, value,  refh="l", refv="b"):
        """
        INPUT:
            refh -- l (left), r (right), c (center)
                    the default is left
            refv -- t (top), b (bottom), B (baseline)
                    the default is bottom.
        """
        if refh == "c": refh=""
        self.append(Text(position, value, refh=refh, refv=refv))

    def ticks(self,
              xmin=None, xmax=None, xscale=1.0, xmult=None, xlen=None, above=True, \
              ymin=None, ymax=None, yscale=1.0, ymult=None, ylen=None, right=True):
        """
        Put a vertical tick mark and label at each position
        on the x-axis between xmin and xmax that when scaled
        is an integer multiple of xmult.  Note that
        xmult must be an integer.  Same with y axis.
        """
        if xmin == None:
            xmin = self.__lower_left[0]
        if xmax == None:
            xmax = self.__upper_right[0]
        if ymin == None:
            ymin = self.__lower_left[1]
        if ymax == None:
            ymax = self.__upper_right[1]
        if xmult== None:
            xmult = int((xmax-xmin)/20.0+1)
        if ymult== None:
            ymult = int((ymax-ymin)/20.0+1)
        if xlen == None:
            xlen = (self.__upper_right[0] - self.__lower_left[0])/60.0
        if ylen == None:
            ylen = (self.__upper_right[1] - self.__lower_left[1])/60.0
        self.append(Ticks(xmin, xmax, xscale, xmult, xlen, above, \
                          ymin, ymax, yscale, ymult, ylen, right))

    def tif(self, density=150):
        if "tif" in self.__done:
            return
        self.ps()
        print "Converting to tif"
        os.system('cd %s; %s -crop 0x0 -density %sx%s "%s.ps" "%s.tif"'%\
                  (tools["GRAPHS_DIR"], tools["CONVERT"], density, density, \
                   self.name(),self.name()))
        print "Done converting"
        self.__done.append("tif")

    def unit(self):
        return self.__unit

    def upper_right(self):
        return self.__upper_right

    def view(self):
        self.ps()
        os.system('cd %s; %s "%s.ps"&'%(tools["GRAPHS_DIR"], tools["PS_VIEWER"],self.name()))






"""
EXAMPLES:

def ex1():
    g=Graph("ex1", (-3,-3),(5,5))
    g.color(0.9,0.9,0.9)
    g.grid()
    g.line_width(0.02)
    g.color(1,0,0)
    g.line((-1,-1),(3,3))
    g.color(1.0,0.5,0.1)
    g.line((-1,-1),(3,0))
    g.color(0,1,0)
    g.line((-2,-1.5),(4,-1.5))
    g.color(0,0,0)
    g.text((0,0), "Hello World")
    points = [Point(x/10.0,(x/10.0)**2) for x in range(-20,20)]
    g.curve(points)
    g.color(1,0,0)
    g.dot((0,0),0.1)
    return g

def ex2():
    g=Graph("ex2", (-3,-5),(3,5))
    g.color(0.7)
    g.grid()
    g.color(0.4)
    g.axes()
    g.color(0,0,0.7)
    g.elliptic_curve([0,0,1,-1,0])
    return g

def tangents():
    g=Graph("tangents",Point(-2,-5),Point(3,5))
    g.color(0.8)
    g.grid()
    g.color(0.3)
    g.axes()
    g.line_width(0.04)
    g.color(0,0,1)
    g.elliptic_curve([0,0,1,-1,0],res=20)
    g.color(1,0,0)
    g.line_width(0.03)
    size = 0.1

    g0 = g.copy("tan0")
    g0.line(Point(-1.5,1.5),Point(2.5,-2.5), "->")
    g0.color(0)
    g0.dot(Point(0,0),size)
    g0.dot(Point(1,-1),size)
    g0.text(Point(2.4,-0.7),"$(1,-1)=t(0,0)$")

    g1 = g.copy("tan1")
    g1.line(Point(0,1),Point(3,-5),"->")
    g1.color(0)
    g1.dot(Point(1,-1),size)
    g1.dot(Point(2,-3),size)
    g1.text(Point(3.6,-2.7), "$(2,-3)=t(1,-1)$")

    g2 = g.copy("tan2")
    delta = 0.5*(Point(2,-3) - Point(21.0/25, -56.0/125))
    g2.line(Point(2,-3)+delta,Point(21.0/25, -56.0/125)-delta, "->")
    g2.color(0)
    g2.dot(Point(21.0/25, -56.0/125), size)
    g2.dot(Point(2,-3),size)
    g2.text(Point(2.6,-0.5), "$\\left(\\frac{21}{25},\\frac{56}{125}\\right) = t(2,-3)$")

    return (g0,g1,g2)
"""
