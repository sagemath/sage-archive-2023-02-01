#!/usr/bin/python
# mplot3d.py
#
# Copyright (C) 2005
# Author:  <jdp@localhost.localdomain>
# Created: 25 Aug 2005
#
# Desc:
#

import random
import math

from matplotlib import lines
from matplotlib import patches
from matplotlib.axes import Axes
from matplotlib import axis
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import normalize
from matplotlib.cm import jet
from matplotlib.transforms import unit_bbox
from matplotlib import cbook
import matplotlib.numerix as nx
import matplotlib.text

import pylab as p
import proj3d

import time

nxa = nx.array
#

def sensible_format_data(self, value):
    if abs(value) > 1e4 or abs(value)<1e-3:
        s = '%1.4e'% value
        return self._formatSciNotation(s)
    else:
        return '%4.3f' % value

def norm_angle(a):
    a = (a+360)%360
    if a > 180: a = a-360
    return a


def text_update_coords(self, renderer):
    """Modified method update_coords from TextWithDash

    I could not understand the original text offset calculations and
    it gave bad results for the angles I was using.  This looks
    better, although the text bounding boxes look a little
    inconsistent
    """

    from matplotlib.numerix import sin, cos, pi, cumsum, dot, asarray, array, \
         where, nonzero, equal, sqrt

    (x, y) = self.get_position()
    dashlength = self.get_dashlength()

    # Shortcircuit this process if we don't have a dash
    if dashlength == 0.0:
        self._mytext.set_position((x, y))
        return

    dashrotation = self.get_dashrotation()
    dashdirection = self.get_dashdirection()
    dashpad = self.get_dashpad()
    dashpush = self.get_dashpush()
    transform = self.get_transform()

    angle = matplotlib.text.get_rotation(dashrotation)

    theta = pi*(angle/180.0+dashdirection-1)
    cos_theta, sin_theta = cos(theta), sin(theta)

    # Compute the dash end points
    # The 'c' prefix is for canvas coordinates
    cxy = array(transform.xy_tup((x, y)))
    cd = array([cos_theta, sin_theta])
    c1 = cxy+dashpush*cd
    c2 = cxy+(dashpush+dashlength)*cd
    (x1, y1) = transform.inverse_xy_tup(tuple(c1))
    (x2, y2) = transform.inverse_xy_tup(tuple(c2))
    self.dashline.set_data((x1, x2), (y1, y2))

    # We now need to extend this vector out to
    # the center of the text area.
    # The basic problem here is that we're "rotating"
    # two separate objects but want it to appear as
    # if they're rotated together.
    # This is made non-trivial because of the
    # interaction between text rotation and alignment -
    # text alignment is based on the bbox after rotation.
    # We reset/force both alignments to 'center'
    # so we can do something relatively reasonable.
    # There's probably a better way to do this by
    # embedding all this in the object's transformations,
    # but I don't grok the transformation stuff
    # well enough yet.
    we = self._mytext.get_window_extent(renderer=renderer)
    w, h = we.width(), we.height()
    off = array([cos_theta*(w/2+2)-1,sin_theta*(h+1)-1])
    off = array([cos_theta*(w/2),sin_theta*(h/2)])
    dir = array([cos_theta,sin_theta])*dashpad
    cw = c2 + off +dir

    self._mytext.set_position(transform.inverse_xy_tup(tuple(cw)))
    # Now set the window extent
    # I'm not at all sure this is the right way to do this.
    we = self._mytext.get_window_extent(renderer=renderer)
    self._window_extent = we.deepcopy()
    self._window_extent.update(((c1[0], c1[1]),), False)

    # Finally, make text align center
    self._mytext.set_horizontalalignment('center')
    self._mytext.set_verticalalignment('center')


class Line3D(lines.Line2D):
    """Make a 2D line pretend to be 3D"""
    def __init__(self, xs,ys,zs, *args, **kwargs):
        lines.Line2D.__init__(self, xs,ys, *args, **kwargs)
        self.xs,self.ys,self.zs = xs,ys,zs

    def draw(self, renderer):
        xs,ys,zs = proj3d.proj_transform(self.xs,self.ys,self.zs, renderer.M)
        self._x,self._y = xs,ys

        lines.Line2D.draw(self, renderer)

class Line3DCollection(LineCollection):
    def __init__(self, segments, *args, **kwargs):
        LineCollection.__init__(self, segments, *args, **kwargs)
        self.segments_3d = segments

    def draw(self, renderer):
        orig_segments = self._segments
        xyslist = [
            proj3d.proj_trans_points(points, renderer.M) for points in
            self.segments_3d]
        segments_2d = [zip(xs,ys) for (xs,ys,zs) in xyslist]
        self._segments = segments_2d
        LineCollection.draw(self, renderer)
        self._segments = orig_segments

class Poly3DCollection(PolyCollection):
    def __init__(self, segments, *args, **kwargs):
        PolyCollection.__init__(self, segments, *args, **kwargs)
        self._zsort = 0
        self.get_vector()

    def get_vector(self):
        """optimise points for projection"""
        si = 0
        ei = 0
        segis = []
        points = []
        for p in self._verts:
            points.extend(p)
            ei = si+len(p)
            segis.append((si,ei))
            si = ei
        xs,ys,zs = zip(*points)
        ones = nx.ones(len(xs))
        self.vec = nx.array([xs,ys,zs,ones])
        self.segis = segis

    def draw(self, renderer):
        #
        #t = time.time()
        txs,tys,tzs = proj3d.proj_transform_vec(self.vec,renderer.M)
        xyslist = [(txs[si:ei],tys[si:ei],tzs[si:ei]) for si,ei in self.segis]
        #print 'simple', time.time()-t

        #t = time.time()
        orig_segments = self._verts
        if 0:
            xyslist = [
                proj3d.proj_trans_points(points, renderer.M) for points in
                self._verts]
        #print 'proj', time.time()-t
        #t = time.time()
        # if required sort by depth (furthest drawn first)
        orig_colors = self._facecolors
        colors = get_colors(orig_colors, len(self._verts))
        #
        if self._zsort:
            z_segments_2d = [(min(zs),zip(xs,ys),c) for (xs,ys,zs),c in zip(xyslist,orig_colors)]
            z_segments_2d.sort()
            z_segments_2d.reverse()
        segments_2d = [s for z,s,c in z_segments_2d]
        colors = [c for z,s,c in z_segments_2d]
        self._verts = segments_2d
        self._facecolors = colors
        #print 'mess', time.time()-t
        #t = time.time()
        PolyCollection.draw(self, renderer)
        #print 'drew', time.time()-t
        self._verts = orig_segments
        self._facecolors = orig_colors

def juggle_axes(xs,ys,zs, dir):
    """Depending on the direction of the plot re-order the axis

    This is so that 2d plots can be plotted along any direction.
    """
    if dir == 'x': return zs,xs,ys
    elif dir == 'y': return xs,zs,ys
    else: return xs,ys,zs

def line_draw(self, renderer):
    """Draw a 2D line as a 3D line"""
    oxs,oys = self.get_xdata(),self.get_ydata()
    xs,ys,zs = juggle_axes(oxs,oys,self.zs,self.dir)
    xs,ys,zs = proj3d.proj_transform(xs,ys,zs, renderer.M)
    self._x = xs
    self._y = ys
    self.old_draw(renderer)
    self._x = oxs
    self._y = oys

def wrap_line(line, zs,dir='z'):
    """Wrap a 2D line so that it draws as a 3D line"""
    line.zs = zs
    line.dir = dir
    line.old_draw = line.draw
    def wrapped_draw(renderer,line=line):
        return line_draw(line,renderer)
    line.draw = wrapped_draw

def set_line_data(line, xs,ys,zs):
    try: line = line[0]
    except: pass
    line.set_data(xs,ys)
    line.zs = zs

def iscolor(c):
    try:
        return (len(c)==4 or len(c)==3) and (type(c[0])==float)
    except (IndexError):
        return None

def get_colors(c, num):
    """Stretch the color argument to provide the required number num"""
    if type(c)==type("string"):
        c = colors.colorConverter.to_rgba(colors)
    if iscolor(c):
        return [c]*num
    elif iscolor(c[0]):
        return c*num
    elif len(c)==num:
        return c[:]
    else:
        raise ValueError, 'unknown color format %s' % c

def zalpha(colors, zs):
    """Modify the alphas of the color list according to depth"""
    colors = get_colors(colors,len(zs))
    norm = normalize(min(zs),max(zs))
    sats = nx.array([1-norm(z)*0.7 for z in zs])
    colors = [(c[0],c[1],c[2],c[3]*s) for c,s in zip(colors,sats)]
    return colors

def patch_draw(self, renderer):
    orig_offsets = self._offsets
    xs,ys = zip(*self._offsets)
    xs,ys,zs = juggle_axes(xs,ys,self.zs,self.dir)
    xs,ys,zs = proj3d.proj_transform(xs,ys,zs, renderer.M)
    # mess with colors
    orig_fcolors = self._facecolors
    orig_ecolors = self._edgecolors
    self._facecolors = zalpha(orig_fcolors,zs)
    self._edgecolors = zalpha(orig_ecolors,zs)

    self._offsets = zip(xs,ys)
    self.old_draw(renderer)
    self._offsets = orig_offsets
    self._facecolors = orig_fcolors
    self._edgecolors = orig_ecolors

def wrap_patch(patch, zs,dir='z',fn=patch_draw):
    patch.zs = zs
    patch.dir = dir
    patch.old_draw = patch.draw
    def wrapped_draw(renderer,patch=patch,fn=fn):
        return fn(patch,renderer)
    patch.draw = wrapped_draw


def draw_linec(self, renderer):
    orig_segments = self._segments
    segments_3d = [[(x,y,z) for (x,y),z in zip(points,zs)]
                   for zs, points in zip(self.zs, self._segments)]
    xyslist = [
        proj3d.proj_trans_points(points, renderer.M) for points in
        segments_3d]
    segments_2d = [zip(xs,ys) for (xs,ys,zs) in xyslist]
    self._segments = segments_2d
    LineCollection.draw(self, renderer)
    self._segments = orig_segments

def draw_polyc(self, renderer):
    orig_segments = self._verts
    # process the list of lists of 2D points held in _verts to generate
    # a list of lists of 3D points
    segments_3d = [[(x,y,z) for (x,y),z in zip(points,self.zs)]
                   for points in self._verts]
    #
    xyslist = [
        proj3d.proj_trans_points(points, renderer.M) for points in
        segments_3d]
    segments_2d = [zip(xs,ys) for (xs,ys,zs) in xyslist]
    self._verts = segments_2d
    PolyCollection.draw(self, renderer)
    self._verts = orig_segments


def text_draw(self, renderer):
    x,y = self.get_position()
    xs,ys,zs = juggle_axes(x,y,self._z,self.dir)
    xs,ys,zs = proj3d.proj_transform(xs,ys,zs, renderer.M)
    self.set_x(xs)
    self.set_y(ys)
    self.old_draw(renderer)
    self.set_x(x)
    self.set_y(y)

def wrap_text(text, zs, dir='z'):
    text._z = zs
    text.dir = dir
    text.old_draw = text.draw
    def wrapped_draw(renderer,text=text):
        return text_draw(text,renderer)
    text.draw = wrapped_draw

def set_text_data(text, x,y,z):
    text._x,text._y,text._z = x,y,z


def draw(text, renderer):
    print 'call draw text', text
    print text.get_visible()
    print 'text "%s"' % text._text
    res = text._get_layout(renderer)
    print res
    text._draw(renderer)

def owrap(text):
    text._draw = text.draw
    def draw_text(renderer,text=text):
        draw(text,renderer)
    text.draw = draw_text

def tick_update_position(tick, x,y,z, angle):
    #
    tick.tick1On = False
    tick.tick2On = False
    tick.tick1line.set_data((x, x),(y,y))
    tick.tick2line.set_data((x, x),(y,y))
    tick.gridline.set_data((x, x),(y,y))
    #
    tick.label1.set_dashlength(8)
    tick.label1.set_dashrotation(angle)
    tick.label1.set_position((x,y))
    tick.label2.set_position((x,y))


class Axis(axis.XAxis):
    def __init__(self, adir, v_intervalx, d_intervalx, axes, *args, **kwargs):
        self.adir = adir
        # data and viewing intervals for this direction
        self.d_interval = d_intervalx
        self.v_interval = v_intervalx
        #
        axis.XAxis.__init__(self, axes, *args, **kwargs)
        self.line = lines.Line2D(xdata=(0,0),ydata=(0,0),
                                 linewidth=0.75,
                                 color=(0,0,0,0),
                                 antialiased=True,
                           )
        #
        # these are the panes which surround the boundary of the view
        self.pane_bg_color = (0.95,0.95,0.95,0.1)
        self.pane_fg_color = (0.9,0.9,0.9,0.5)
        self.pane = patches.Polygon([],
                                    alpha=0.2,
                                    facecolor=self.pane_fg_color,
                                    edgecolor=self.pane_fg_color)
        #
        self.axes._set_artist_props(self.line)
        self.axes._set_artist_props(self.pane)
        self.gridlines = Line3DCollection([])
        self.axes._set_artist_props(self.gridlines)
        self.axes._set_artist_props(self.label)
        self.label._transform = self.axes.transData


    def get_tick_positions(self):
        majorTicks = self.get_major_ticks()
        majorLocs = self.major.locator()
        self.major.formatter.set_locs(majorLocs)
        majorLabels = [self.major.formatter(val, i) for i, val in enumerate(majorLocs)]
        return majorLabels,majorLocs

    def get_major_ticks(self):
        ticks = axis.XAxis.get_major_ticks(self)
        for t in ticks:
            def update_coords(renderer,self=t.label1):
                return text_update_coords(self, renderer)
            # Text overrides setattr
            t.label1.__dict__['update_coords'] = update_coords
            t.tick1line.set_transform(self.axes.transData)
            t.tick2line.set_transform(self.axes.transData)
            t.gridline.set_transform(self.axes.transData)
            t.label1.set_transform(self.axes.transData)
            t.label2.set_transform(self.axes.transData)
        return ticks

    def set_pane_fg(self, xys):
        self.pane.xy = xys
        self.pane.set_edgecolor(self.pane_fg_color)
        self.pane.set_facecolor(self.pane_fg_color)
        self.pane.set_alpha(self.pane_fg_color[-1])

    def set_pane_bg(self, xys):
        self.pane.xy = xys
        self.pane.set_edgecolor(self.pane_bg_color)
        self.pane.set_facecolor(self.pane_bg_color)
        self.pane.set_alpha(self.pane_bg_color[-1])

    def draw(self, renderer):
        #
        renderer.open_group('axis3d')
        ticklabelBoxes = []
        ticklabelBoxes2 = []

        # code from XAxis
        majorTicks = self.get_major_ticks()
        majorLocs = self.major.locator()
        self.major.formatter.set_locs(majorLocs)
        majorLabels = [self.major.formatter(val, i)
                       for i, val in enumerate(majorLocs)]
        #
        minx,maxx,miny,maxy,minz,maxz = self.axes.get_w_lims()

        interval = self.get_view_interval()
        # filter locations here so that no extra grid lines are drawn
        majorLocs = [loc for loc in majorLocs if interval.contains(loc)]
        # these will generate spacing for labels and ticks
        dx = (maxx-minx)/12
        dy = (maxy-miny)/12
        dz = (maxz-minz)/12

        # stretch the boundary slightly so that the ticks have a better fit
        minx,maxx,miny,maxy,minz,maxz = (
            minx-dx/4,maxx+dx/4,miny-dy/4,maxy+dy/4,minz-dz/4,maxz+dz/4)

        # generate the unit_cubes and transformed unit_cubes from the stretched
        # limits
        vals = minx,maxx,miny,maxy,minz,maxz
        uc = self.axes.unit_cube(vals)
        tc = self.axes.tunit_cube(vals,renderer.M)
        #
        # these are flags which decide whether the axis should be drawn
        # on the high side (ie on the high side of the paired axis)
        xhigh = tc[1][2]>tc[2][2]
        yhigh = tc[3][2]>tc[2][2]
        zhigh = tc[0][2]>tc[2][2]
        #
        aoff = 0

        #
        if self.adir == 'x':
            lx = (minx+maxx)/2
            if xhigh:
                # xaxis at front
                self.set_pane_fg([tc[0],tc[1],tc[5],tc[4]])
                to = tc[3]
                te = tc[2]
                xyz = [(x,maxy,minz) for x in majorLocs]
                nxyz = [(x,miny,minz) for x in majorLocs]
                lxyz = [(x,miny,maxz) for x in majorLocs]
                aoff = -90

                ly = maxy + dy
                lz = minz - dz
            else:
                self.set_pane_bg([tc[3],tc[2],tc[6],tc[7]])
                to = tc[0]
                te = tc[1]
                xyz = [(x,miny,minz) for x in majorLocs]
                nxyz = [(x,maxy,minz) for x in majorLocs]
                lxyz = [(x,maxy,maxz) for x in majorLocs]
                aoff = 90

                ly = miny - dy
                lz = minz - dz
        elif self.adir == 'y':
            # cube 3 is minx,maxy,minz
            # cube 2 is maxx,maxy,minz
            ly = (maxy+miny)/2
            if yhigh:
                # yaxis at front
                self.set_pane_fg([tc[0],tc[3],tc[7],tc[4]])
                to = tc[1]
                te = tc[2]
                xyz = [(maxx,y,minz) for y in majorLocs]
                nxyz = [(minx,y,minz) for y in majorLocs]
                lxyz = [(minx,y,maxz) for y in majorLocs]
                aoff = 90

                #
                lx = maxx + dx
                lz = minz - dz

            else:
                # yaxis at back
                self.set_pane_bg([tc[1],tc[5],tc[6],tc[2]])
                to = tc[0]
                te = tc[3]
                xyz = [(minx,y,minz) for y in majorLocs]
                nxyz = [(maxx,y,minz) for y in majorLocs]
                lxyz = [(maxx,y,maxz) for y in majorLocs]
                aoff = -90
                #
                lx = minx - dx
                lz = minz - dz

        elif self.adir == 'z':
            nxyz = None
            self.set_pane_bg([tc[0],tc[1],tc[2],tc[3]])
            aoff = -90
            lz = (maxz+minz)/2
            if xhigh and yhigh:
                to = tc[1]
                te = tc[5]
                xyz = [(maxx,miny,z) for z in majorLocs]
                nxyz = [(minx,miny,z) for z in majorLocs]
                lxyz = [(minx,maxy,z) for z in majorLocs]
                #
                lx = maxx + dx
                ly = miny - dy
            elif xhigh and not yhigh:
                to = tc[2]
                te = tc[6]
                xyz = [(maxx,maxy,z) for z in majorLocs]
                nxyz = [(maxx,miny,z) for z in majorLocs]
                lxyz = [(minx,miny,z) for z in majorLocs]

                lx = maxx + dx
                ly = maxy + dy
            elif yhigh and not xhigh:
                to = tc[0]
                te = tc[4]
                xyz = [(minx,miny,z) for z in majorLocs]
                nxyz = [(minx,maxy,z) for z in majorLocs]
                lxyz = [(maxx,maxy,z) for z in majorLocs]
                lx = minx - dx
                ly = miny - dy
            else:
                to = tc[3]
                te = tc[7]
                xyz = [(minx,maxy,z) for z in majorLocs]
                nxyz = [(maxx,maxy,z) for z in majorLocs]
                lxyz = [(maxx,miny,z) for z in majorLocs]
                lx = minx - dx
                ly = maxy + dy

        #
        tlx,tly,tlz = proj3d.proj_transform(lx,ly,lz, renderer.M)
        self.label.set_position((tlx,tly))
        #print self.label._text, lx,ly, tlx,tly
        #
        # SAGE: don't want
        #self.pane.draw(renderer)
        #
        self.line.set_transform(self.axes.transData)
        self.line.set_data((to[0],te[0]),(to[1],te[1]))
        self.line.draw(renderer)

        angle = norm_angle(math.degrees(math.atan2(te[1]-to[1],te[0]-to[0])))
        #
        # should be some other enabler here...
        if len(self.label._text)>1:
            if abs(angle)>90 and self.adir != 'z':
                self.label.set_rotation(angle+180)
            else:
                self.label.set_rotation(angle)
        self.label.draw(renderer)
        #
        angle = angle + aoff

        points = proj3d.proj_points(xyz,renderer.M)
        if nxyz:
            tnxyz = proj3d.proj_points(nxyz,renderer.M)
            tlxyz = proj3d.proj_points(lxyz,renderer.M)
            lines = zip(xyz,nxyz,lxyz)
            self.gridlines.segments_3d = lines
            #self.gridlines._colors = [(0.9,0.9,0.9,1)]*len(lines)
            self.gridlines._colors = [(0.92,0.92,0.92,1.0)]*len(lines)
            self.gridlines.draw(renderer)

        seen = {}
        interval = self.get_view_interval()
        # SAGE - commented out since suddenly broken.
##         for tick, loc, (x,y,z), label in zip(majorTicks, majorLocs, points, majorLabels):
##             if tick is None: continue
##             if not interval.contains(loc): continue
##             seen[loc] = 1
##             tick_update_position(tick, x,y,z, angle=angle)
##             tick.set_label1(label)
##             tick.set_label2(label)
##             tick.draw(renderer)
##             if tick.label1On:
##                 extent = tick.label1.get_window_extent(renderer)
##                 ticklabelBoxes.append(extent)
##             if tick.label2On:
##                 extent = tick.label2.get_window_extent(renderer)
##                 ticklabelBoxes2.append(extent)
        #
        renderer.close_group('axis3d')


    def get_view_interval(self):
        """return the Interval instance for this axis view limits
        """
        return self.v_interval()

    def get_data_interval(self):
        'return the Interval instance for this axis data limits'
        return self.d_interval()

class Axes3D(Axes):
    def __init__(self, fig, rect=[0.05, 0.05, 0.9, 0.9],
                 elev=30, azim=-60, *args, **kwargs):
        #
        figmanager = p.get_current_fig_manager()
        #self.toolbar = figmanager.toolbar

        self.xy_viewLim = unit_bbox()
        self.zz_viewLim = unit_bbox()
        self.xy_dataLim = unit_bbox()
        self.zz_dataLim = unit_bbox()
        # inihibit autoscale_view until the axises are defined
        # they can't be defined until Axes.__init__ has been called
        self._ready = 0
        Axes.__init__(self, fig, rect,
                      frameon=True,
                      xticks=[], yticks=[], *args, **kwargs)
        #
        self._ready = 1
        self.view_init(elev, -azim)
        self.create_axes()
        self.set_top_view()

    def set_top_view(self):
        # this happens to be the right view for the viewing coordinates
        # moved up and to the left slightly to fit labels and axes
        dwl = 0.95/self.dist
        dw = 0.9/self.dist
        self.set_xlim(-dwl,dw)
        self.set_ylim(-dwl,dw)

    def create_axes(self):
        self.w_xaxis = Axis('x',self.xy_viewLim.intervalx,
                            self.xy_dataLim.intervalx, self)
        self.w_yaxis = Axis('y',self.xy_viewLim.intervaly,
                            self.xy_dataLim.intervaly, self)
        self.w_zaxis = Axis('z',self.zz_viewLim.intervalx,
                            self.zz_dataLim.intervalx, self)


    def unit_cube(self,vals=None):
        minx,maxx,miny,maxy,minz,maxz = vals or self.get_w_lims()
        xs,ys,zs = ([minx,maxx,maxx,minx,minx,maxx,maxx,minx],
                    [miny,miny,maxy,maxy,miny,miny,maxy,maxy],
                    [minz,minz,minz,minz,maxz,maxz,maxz,maxz])
        return zip(xs,ys,zs)

    def tunit_cube(self,vals=None,M=None):
        if M is None: M = self.M
        xyzs = self.unit_cube(vals)
        tcube = proj3d.proj_points(xyzs,M)
        return tcube

    def tunit_edges(self, vals=None,M=None):
        tc = self.tunit_cube(vals,M)
        edges = [(tc[0],tc[1]),
                 (tc[1],tc[2]),
                 (tc[2],tc[3]),
                 (tc[3],tc[0]),

                 (tc[0],tc[4]),
                 (tc[1],tc[5]),
                 (tc[2],tc[6]),
                 (tc[3],tc[7]),

                 (tc[4],tc[5]),
                 (tc[5],tc[6]),
                 (tc[6],tc[7]),
                 (tc[7],tc[4])]
        return edges

    def draw(self, renderer):
        #print 'elev', self.elev, 'azim', self.azim

        # draw the background patch
        self.axesPatch.draw(renderer)
        self._frameon = False

        # add the projection matrix to the renderer
        self.M = self.get_proj()
        renderer.M = self.M
        self.set_top_view()
        self.w_xaxis.draw(renderer)
        self.w_yaxis.draw(renderer)
        self.w_zaxis.draw(renderer)
        Axes.draw(self, renderer)

    def update_datalim(self, xys):
        pass

    def update_datalim_numerix(self, x, y):
        pass

    def auto_scale_xyz(self, X,Y,Z=None,had_data=None):
        return
        x,y,z = X,Y,Z
        try:
            x,y = X.flat,Y.flat
            if Z is not None:
                z = Z.flat
        except AttributeError:
            pass

        self.xy_dataLim.update_numerix(x, y, not had_data)
        if z is not None:
            self.zz_dataLim.update_numerix(z, z, not had_data)
        self.autoscale_view()

    def autoscale_view(self):
        return
        self.set_top_view()
        if not self._ready: return

        if not self._autoscaleon: return
        locator = self.w_xaxis.get_major_locator()
        #print 'auto', locator.autoscale()
        self.set_w_xlim(locator.autoscale())
        locator = self.w_yaxis.get_major_locator()
        self.set_w_ylim(locator.autoscale())
        locator = self.w_zaxis.get_major_locator()
        self.set_w_zlim(locator.autoscale())

    def get_w_lims(self):
        minx,maxx = self.get_w_xlim()
        miny,maxy = self.get_w_ylim()
        minz,maxz = self.get_w_zlim()
        return minx,maxx,miny,maxy,minz,maxz

    def set_w_zlim(self, *args, **kwargs):
        gl,self.get_xlim = self.get_xlim,self.get_w_zlim
        vl,self.viewLim = self.viewLim,self.zz_viewLim
        vmin,vmax = Axes.set_xlim(self, *args, **kwargs)
        self.get_xlim = gl
        self.viewLim = vl
        return vmin,vmax

    def set_w_xlim(self, *args, **kwargs):
        gl,self.get_xlim = self.get_xlim,self.get_w_xlim
        vl,self.viewLim = self.viewLim,self.xy_viewLim
        vmin,vmax = Axes.set_xlim(self, *args, **kwargs)
        self.get_xlim = gl
        self.viewLim = vl
        return vmin,vmax

    def set_w_ylim(self, *args, **kwargs):
        gl,self.get_ylim = self.get_ylim,self.get_w_ylim
        vl,self.viewLim = self.viewLim,self.xy_viewLim
        vmin,vmax = Axes.set_ylim(self, *args, **kwargs)
        self.get_ylim = gl
        self.viewLim = vl
        return vmin,vmax

    def get_w_zlim(self):
        return self.zz_viewLim.intervalx().get_bounds()

    def get_w_xlim(self):
        return self.xy_viewLim.intervalx().get_bounds()

    def get_w_ylim(self):
        return self.xy_viewLim.intervaly().get_bounds()


    def get_affine(self, z):
        a,b,zfx,ttx = self.M[0]
        c,d,zfy,tty = self.M[1]
        tx = z*zfx + ttx
        ty = z*zfy + tty
        return affine(a,b,c,d,tx,ty)

    def get_proj(self):

        relev,razim = nx.pi * self.elev/180, nx.pi * self.azim/180

        xmin,xmax = self.get_w_xlim()
        ymin,ymax = self.get_w_ylim()
        zmin,zmax = self.get_w_zlim()

        # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
        worldM = proj3d.world_transformation(xmin,xmax,
                                             ymin,ymax,
                                             zmin,zmax)

        # look into the middle of the new coordinates
        R = nxa([0.5,0.5,0.5])
        #
        xp = R[0] + nx.cos(razim)*nx.cos(relev)*self.dist
        yp = R[1] + nx.sin(razim)*nx.cos(relev)*self.dist
        zp = R[2] + nx.sin(relev)*self.dist

        E = nxa((xp, yp, zp))
        if abs(relev) > nx.pi/2:
            # upside down
            V = nxa((0,0,-1))
        else:
            V = nxa((0,0,1))
        zfront,zback = -self.dist,self.dist

        viewM = proj3d.view_transformation(E,R,V)
        perspM = proj3d.persp_transformation(zfront,zback)
        M0 = nx.matrixmultiply(viewM,worldM)
        M = nx.matrixmultiply(perspM,M0)
        return M

    def default_coord(self, kwargs):
        v = 'z'
        zval = 0
        if kwargs.has_key('zval'):
            zval = kwargs['zval']
            del kwargs['zval']
            v = 'z'
        if kwargs.has_key('yval'):
            zval = kwargs['yval']
            del kwargs['yval']
            v = 'y'
        if kwargs.has_key('xval'):
            zval = kwargs['xval']
            del kwargs['xval']
            v = 'x'
        return zval,v

    def plot(self, *args, **kwargs):
        had_data = self.has_data()

        zval,zdir = self.default_coord(kwargs)

        lines = Axes.plot(self, *args, **kwargs)
        line = lines[0]
        xs = line.get_xdata()
        ys = line.get_ydata()
        zs = nx.array([zval for x in xs])
        wrap_line(line,zs,dir=zdir)
        xs,ys,zs = juggle_axes(xs,ys,zs,zdir)
        #
        self.auto_scale_xyz(xs,ys,zs, had_data)


    def plot3D(self, xs, ys, zs, *args, **kwargs):
        had_data = self.has_data()
        lines = Axes.plot(self, xs,ys, *args, **kwargs)
        if len(lines)==1:
            line = lines[0]
            wrap_line(line, zs)
        #
        self.auto_scale_xyz(xs,ys,zs, had_data)
        return lines

    def plot_surface(self, X, Y, Z, *args, **kwargs):
        had_data = self.has_data()

        rows,cols = Z.shape
        tX,tY,tZ = p.transpose(X),p.transpose(Y),p.transpose(Z)
        maxr = rows
        div = cbook.popd(kwargs, 'div', 20)
        # this should be smarter...
        rstride = rows/div
        maxc = cols
        cstride = cols/div
        polys = []
        boxes = []
        for rs in p.arange(0,maxr,rstride):
            for cs in p.arange(0,maxc,cstride):
                ps = []
                corners = []
                for a,ta in [(X,tX),(Y,tY),(Z,tZ)]:
                    ztop = a[rs][cs:min(maxc-1,cs+cstride)]
                    zleft = ta[min(maxc-1,cs+cstride)][rs:min(maxr-1,rs+rstride)]
                    zbase = a[min(maxr-1,rs+rstride)][cs:min(maxc-1,cs+cstride):]
                    zbase = zbase[::-1]
                    zright = ta[cs][rs:min(maxr-1,rs+rstride):]
                    zright = zright[::-1]
                    corners.append([ztop[0],ztop[-1],zbase[0],zbase[-1]])
                    z = p.concatenate((ztop,zleft,zbase,zright))
                    ps.append(z)
                boxes.append(map(nxa,zip(*corners)))
                polys.append(zip(*ps))
        #

        #
        #
        lines = []
        shade = []
        for box in boxes:
            n = proj3d.cross(box[0]-box[1],
                         box[0]-box[2])
            n = n/proj3d.mod(n)*5
            shade.append(proj3d.dot(n,[-1,-1,0]))
            lines.append((box[0],n+box[0]))
        #
        #self.add_lines(lines, colors=[(1.0,0,0,1.0)]*len(lines))
        #
        color = nxa([0,0,1,1])
        norm = normalize(min(shade+[0]),max(shade+[0]))
        if 0:
            cmap = jet
            colors = [cmap(norm(v)) for v in shade]
        else:
            colors = [color * (0.5+norm(v)*0.5) for v in shade]
            for c in colors: c[3] = 1
            #colors = [(c[:3]+[1]) for c in colors]
        if 0:
            polyc = Poly3DCollection(boxes, facecolors=colors, *args, **kwargs)
            polyc._zsort = 1
            self.add_collection(polyc)
        else:
            polyc = Poly3DCollection(polys, facecolors=colors, *args, **kwargs)
            polyc._zsort = 1
            self.add_collection(polyc)

        self.auto_scale_xyz(X,Y,Z, had_data)
        return polyc

    def plot_wireframe(self, X, Y, Z, *args, **kwargs):
        nrow = cbook.popd(kwargs, "nrow", 1)
        ncol = cbook.popd(kwargs, "ncol", 1)

        had_data = self.has_data()
        rows,cols = Z.shape

        tX,tY,tZ = p.transpose(X),p.transpose(Y),p.transpose(Z)

        rii = [i for i in range(0,rows,nrow)]+[rows-1]
        cii = [i for i in range(0,cols,ncol)]+[cols-1]
        xlines = [X[i] for i in rii]
        ylines = [Y[i] for i in rii]
        zlines = [Z[i] for i in rii]
        #
        txlines = [tX[i] for i in cii]
        tylines = [tY[i] for i in cii]
        tzlines = [tZ[i] for i in cii]
        #
        lines = [zip(xl,yl,zl) for xl,yl,zl in zip(xlines,ylines,zlines)]
        lines += [zip(xl,yl,zl) for xl,yl,zl in zip(txlines,tylines,tzlines)]
        linec = self.add_lines(lines, *args, **kwargs)

        self.auto_scale_xyz(X,Y,Z, had_data)
        return linec

    def contour3D(self, X, Y, Z, *args, **kwargs):
        had_data = self.has_data()
        #levels, colls = self.contour(X, Y, Z, *args, **kwargs)
        C = self.contour(X, Y, Z, *args, **kwargs)
        levels, colls = (C.levels, C.collections)
        for z,linec in zip(levels,colls):
            zl = []
            for s in linec._segments:
                zl.append([z] * len(s))
            wrap_patch(linec, zl, fn=draw_linec)
        self.auto_scale_xyz(X,Y,Z, had_data)
        return levels,colls

    def contourf3D(self, X, Y, Z, *args, **kwargs):
        had_data = self.has_data()

        #levels, colls = self.contourf(X, Y, Z, 20)
        C = self.contourf(X, Y, Z, *args, **kwargs)
        levels, colls = (C.levels, C.collections)
        print len(levels),len(colls)
        for z1,z2,linec in zip(levels,levels[1:],colls):
            zs = [z1] * (len(linec._verts[0])/2)
            zs += [z2] * (len(linec._verts[0])/2)
            wrap_patch(linec, zs, fn=draw_polyc)
        self.auto_scale_xyz(X,Y,Z, had_data)
        return levels,colls


    def scatter3D(self, xs, ys, zs, *args, **kwargs):
        had_data = self.has_data()
        patches = Axes.scatter(self,xs,ys,*args,**kwargs)
        wrap_patch(patches, zs)
        #
        self.auto_scale_xyz(xs,ys,zs, had_data)

    def add_lines(self, lines, *args, **kwargs):
        linec = Line3DCollection(lines, *args, **kwargs)
        self.add_collection(linec)

    def text3D(self, x,y,z,s, *args, **kwargs):
        text = Axes.text(self,x,y,s,*args,**kwargs)
        wrap_text(text,z)
        return text

    def view_init(self, elev, azim):
        self.dist = 10
        self.elev = elev
        self.azim = azim
        #
        self.button_pressed = None
        #self.figure.canvas.mpl_connect('motion_notify_event', self.on_move)

        #self.figure.canvas.mpl_connect('button_press_event', self.button_press)
        #self.figure.canvas.mpl_connect('button_release_event', self.button_release)

    def button_press(self, event):
        if event.button == 2:
            self.button_pressed = 2
            self.sx,self.sy = event.xdata,event.ydata

    def button_release(self, event):
        if event.button == 2:
            self.button_pressed = None
            #self.figure.canvas.mpl_disconnect(self.bind_id)



    def format_xdata(self, x):
        """
        Return x string formatted.  This function will use the attribute
        self.fmt_xdata if it is callable, else will fall back on the xaxis
        major formatter
        """
        try: return self.fmt_xdata(x)
        except TypeError:
            fmt = self.w_xaxis.get_major_formatter()
            return sensible_format_data(fmt,x)

    def format_ydata(self, y):
        """
        Return y string formatted.  This function will use the attribute
        self.fmt_ydata if it is callable, else will fall back on the yaxis
        major formatter
        """
        try: return self.fmt_ydata(y)
        except TypeError:
            fmt = self.w_yaxis.get_major_formatter()
            return sensible_format_data(fmt,y)

    def format_zdata(self, z):
        """
        Return y string formatted.  This function will use the attribute
        self.fmt_ydata if it is callable, else will fall back on the yaxis
        major formatter
        """
        try: return self.fmt_zdata(z)
        except (AttributeError,TypeError):
            fmt = self.w_zaxis.get_major_formatter()
            return sensible_format_data(fmt,z)


    def format_coord(self, xd, yd):
        """Given the 2D view coordinates attempt to guess a 3D coordinate

        Looks for the nearest edge to the point and then assumes that the point is
        at the same z location as the nearest point on the edge.
        """
        p = (xd,yd)
        edges = self.tunit_edges()
        #lines = [proj3d.line2d(p0,p1) for (p0,p1) in edges]
        ldists = [(proj3d.line2d_seg_dist(p0,p1,p),i) for i,(p0,p1) in enumerate(edges)]
        ldists.sort()
        # nearest edge
        edgei = ldists[0][1]
        #
        p0,p1 = edges[edgei]

        # scale the z value to match
        x0,y0,z0 = p0
        x1,y1,z1 = p1
        d0 = nx.hypot(x0-xd,y0-yd)
        d1 = nx.hypot(x1-xd,y1-yd)
        dt = d0+d1
        z = d1/dt * z0 + d0/dt * z1
        #print 'mid', edgei, d0, d1, z0, z1, z

        x,y,z = proj3d.inv_transform(xd,yd,z,self.M)

        xs = self.format_xdata(x)
        ys = self.format_ydata(y)
        zs = self.format_ydata(z)
        return  'x=%s, y=%s z=%s'%(xs,ys,zs)

    def on_move(self, event):
        if event.inaxes != self:
            return
        #
        if not self.button_pressed:
            s = event.inaxes.format_coord(event.xdata, event.ydata)
            #self.toolbar.set_message(s)
            return
        #

        # get the x and y pixel coords

        x, y = event.xdata, event.ydata
        dx,dy = x-self.sx,y-self.sy
        self.sx,self.sy = x,y
        if dx == 0 and dy == 0: return

        #
        x0,x1 = self.get_xlim()
        y0,y1 = self.get_ylim()
        w = (x1-x0)
        h = (y1-y0)
        #
        self.elev = norm_angle(self.elev - (dy/h)*180)
        self.azim = norm_angle(self.azim - (dx/w)*180)
        #self.toolbar.set_message('azimuth=%d deg, elevation=%d deg ' % (self.azim, self.elev))
        self.get_proj()
        self.figure.canvas.draw()


    def set_xlabel(self, xlabel, fontdict=None, **kwargs):
        #par = cbook.popd(kwargs, 'par',None)
        #label.set_par(par)
        #
        label = self.w_xaxis.get_label()
        label.set_text(xlabel)
        if fontdict is not None: label.update(fontdict)
        label.update(kwargs)
        return label

    def set_ylabel(self, ylabel, fontdict=None, **kwargs):
        label = self.w_yaxis.get_label()
        label.set_text(ylabel)
        if fontdict is not None: label.update(fontdict)
        label.update(kwargs)
        return label

    def set_zlabel(self, zlabel, fontdict=None, **kwargs):
        label = self.w_zaxis.get_label()
        label.set_text(zlabel)
        if fontdict is not None: label.update(fontdict)
        label.update(kwargs)
        return label


def test1():
    fig = p.figure()
    #x,y,z = nx.mlab.rand(3,200)
    ax = Axes3D(fig)
    #
    #
    n = 200
    for c,zl,zh in [('r',-50,-25),('b',-30,-5)]:
        xs,ys,zs = zip(*
                       [(random.randrange(23,32),
                         random.randrange(100),
                         random.randrange(zl,zh)
                         ) for i in range(n)])
        ax.scatter3D(xs,ys,zs, c=c)
    #
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    return ax

    fig.add_axes(ax)
    return fig
    p.show()


def test2():
    fig = p.figure()
    ax = Axes3D(fig)
    #
    xs = nx.arange(0,20*p.pi,0.1)
    ys = 50*nx.sin(xs/10)
    ax.plot(xs,ys,zval=0)
    ax.plot(xs,ys,xval=0)
    ax.plot(xs,ys,yval=0)
    ax.plot(xs,ys,xval=max(xs))
    ax.plot(xs,ys,yval=max(ys))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #
    fig.add_axes(ax)
    p.show()

def test3():
    fig = p.figure()
    ax = Axes3D(fig)

    Cos,Sin,E,Pi = p.cos,p.sin,p.exp(1),p.pi
    xl,yl,zl = [],[],[]
    for u in p.arange(0,15,0.1):
        xs,ys,zs = [],[],[]
        for v in p.arange(0,2*p.pi+0.2,0.2):
            x = 2*(1 - E**(u/(6*Pi)))*Cos(u)*Cos(v/2)**2
            y = 2*(-1 + E**(u/(6*Pi)))*Cos(v/2)**2*Sin(u)
            z = 1 - E**(u/(3*Pi)) - Sin(v) + E**(u/(6*Pi))*Sin(v)
            xs.append(x)
            ys.append(y)
            zs.append(z)
        xl.append(xs)
        yl.append(ys)
        zl.append(zs)
    X,Y,Z = p.array(xl),p.array(yl),p.array(zl)
    print X.shape,Y.shape,Z.shape
    #ax.plot_wireframe(p.array(xl),p.array(yl),p.array(zl),ncol=1,nrow=1)


    if 1:
        ax.plot3D(p.array(xs),p.array(ys),p.array(zs))
        xs,ys,zs = [],[],[]
        for v in p.arange(0,2*p.pi+0.1,0.1):
            for u in p.arange(0,15,0.1):
                x = 2*(1 - E**(u/(6*Pi)))*Cos(u)*Cos(v/2)**2
                y = 2*(-1 + E**(u/(6*Pi)))*Cos(v/2)**2*Sin(u)
                z = 1 - E**(u/(3*Pi)) - Sin(v) + E**(u/(6*Pi))*Sin(v)
                xs.append(x)
                ys.append(y)
                zs.append(z)
        ax.plot3D(p.array(xs),p.array(ys),p.array(zs))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.add_axes(ax)
    p.savefig("surface2.png")
    return fig

if __name__ == "__main__":
    test3()





