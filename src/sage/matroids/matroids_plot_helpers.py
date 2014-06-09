r"""
Helper functions for plotting the geometric representation of matroids


AUTHORS:

- Jayant Apte (2014-06-06): initial version


EXAMPLES::

"""

#*****************************************************************************
#       Copyright (C) 2013 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import scipy
import scipy.interpolate
import numpy as np
import operator
from sage.structure.sage_object import SageObject
from sage.plot.all import Graphics, line, text, polygon2d, point,points
from sage.sets.set import Set

def initial_triangle(M,B1,nB1,lps):
    """
    Return points on and off the triangle and lines to be drawn for a rank 3 matroid
    
    INPUT:
    
    - ``M`` -- A matroid.
    - ``B1``-- A list of labels of groundset elements of M that corresponds to a basis of matroid
    returned by ``M.simplify()``.
    - ``nB1``-- A list of labes of elements in the ground set of M that corresponds to ``M.simplify.groundset() \ B1``. 
    - ``lps``-- A list of labels of elements in the ground set of matroid M that are loops.
    
    OUTPUT:
    
    A tuple containing 4 elements in this order:
    
    1. A dictionary containing 2-tuple (x,y) co-ordinates with ``M.simplify.groundset()`` elements that can be placed on the sides of the triangle as keys
    2. A list of 3 lists of elements of ``M.simplify.groundset()`` that can be placed on the 3 sides of the triangle
    3. A list of elements of `M.simplify.groundset()`` that cane be placed inside the triangle in the geometric representation
    4. A list of lists of elements of ``M.simplify.groundset()`` that correspond to lines in the geometric representation other than the sides of the triangle
    
    EXAMPLES::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: M=Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0],[0, 1, 0, 1, 0, 1, 1,0],[0, 0, 1, 1, 1, 0, 1,0]])
        sage: N=M.simplify()
        sage: B1=list(N.basis())
        sage: nB1=list(set(M.simplify().groundset())-set(B1))
        sage: pts,trilines,nontripts,curvedlines=matroids_plot_helpers.initial_triangle(M,B1,nB1,M.loops())
        sage: print pts
        {1: (1.0, 0.0), 2: (1.5, 1.0), 3: (0.5, 1.0), 4: (0, 0), 5: (1, 2), 6: (2, 0)}
        sage: print trilines
        [[3, 4, 5], [2, 5, 6], [1, 4, 6]]
        sage: print nontripts
        [0]
        sage: print curvedlines
        [[0, 1, 5], [0, 2, 4], [0, 3, 6], [1, 2, 3], [1, 4, 6], [2, 5, 6], [3, 4, 5]]

    
    .. NOTE::

            This method does NOT do any checks. 
            
    """ 
    tripts = [(0, 0),(1, 2),(2, 0)]
    pts = {}
    j=0
    for i in B1:
        pts[i] = tripts[j]
        j = j + 1
    #nB1=[j for j in gnd if j not in B1]
    pairs=[[0, 1],[1,2],[0,2]]
    L1 = []
    L2 = []
    L3 = []
    for i in nB1:
        if M.is_dependent([i,B1[pairs[0][0]],B1[pairs[0][1]]]):
            # Add to L1 
            L1.append(i)
        elif M.is_dependent([i,B1[pairs[1][0]],B1[pairs[1][1]]]):
            # Add to L2  
            L2.append(i)
        elif M.is_dependent([i,B1[pairs[2][0]],B1[pairs[2][1]]]):
            # Add to L3 
            L3.append(i)
    L = [L1,L2,L3]
    lines = [] #the list of lines
    for i in range(1,len(L)+1):
        lines.append([B1[pairs[i-1][0]]])
        lines[i-1].extend(L[i-1])
        lines[i-1].extend([B1[pairs[i-1][1]]])
    # place triangle and L1,L2,L3
    for i in L:# loop over megalist
        interval = 1/float(len(i)+1)
        pt1 = list(tripts[pairs[L.index(i)][0]])
        pt2 = list(tripts[pairs[L.index(i)][1]])
        for j in range(1,len(i)+1): 
            #Loop over L1,L2,L3
            cc = interval*j
            pts[i[j-1]] = (cc*pt1[0]+(1-cc)*pt2[0],cc*pt1[1]+(1-cc)*pt2[1])
    trilines = [list(set(x)) for x in lines if len(x)>=3]
    curvedlines = [list(set(list(x)).difference(set(lps))) for x in M.flats(2) if set(list(x)) not in trilines if len(list(x))>=3]
    nontripts = [i for i in nB1 if i not in pts.keys()] 
    return pts,trilines,nontripts,curvedlines

def trigrid(tripts):
    """
    Return a grid of 4 points inside given 3 points as a list
    
   
    INPUT:
    
    - ``tripts`` -- A list of 3 lists of the form [x,y] where x and y are the cartesian co-ordinates of a point
    
    OUTPUT:
    
    A list of lists containing 4 points in following order:
    
    1. Barycenter of 3 input points 
    2,3,4. Barycenters of 1. with 3 different 2-subsets of input points respectively  
    
    EXAMPLES::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: points=matroids_plot_helpers.trigrid([[2,1],[4,5],[5,2]])
        sage: print points
        [[3.6666666666666665, 2.6666666666666665], [3.222222222222222, 2.888888888888889], [4.222222222222222, 3.222222222222222], [3.5555555555555554, 1.8888888888888886]]
    
    .. NOTE::

            This method does NOT do any checks.
            
    """ 
    n = 0
    pairs = [[0,1],[1,2],[0,2]]
    cpt = list((float(tripts[0][0]+tripts[1][0]+tripts[2][0])/3,float(tripts[0][1]+tripts[1][1]+tripts[2][1])/3))
    grid = [cpt]
    for p in pairs:
        pt = list((float(tripts[p[0]][0]+tripts[p[1]][0]+cpt[0])/3,float(tripts[p[0]][1]+tripts[p[1]][1]+cpt[1])/3))
        grid.append(pt)
    return grid
    
def addnontripts(tripts_labels,nontripts_labels,ptsdict):
    """
    Return modified ``ptsdict`` with additional keys and values corresponding to ``nontripts``
    
   
    INPUT:
    
    - ``tripts`` -- A list of 3 labels that are to be placed on vertices of the triangle
    - ``ptsdict`` -- A dictionary (at least) containing labels in ``tripts`` as keys and their (x,y) position as values
    - ``nontripts``-- A list of labels whose corresponding points are to be placed inside the triangle
    
    OUTPUT:
    
    A dictionary containing labels in ``tripts`` as keys and their (x,y) position as values allong with all keys and respective values in ``ptsdict`` 
    
    EXAMPLES::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: ptsdict={'a':(0,0),'b':(1,2),'c':(2,0)}
        sage: ptsdict_1=matroids_plot_helpers.addnontripts(['a','b','c'],['d','e','f'],ptsdict)
        sage: print ptsdict_1
        {'a': (0, 0), 'c': (2, 0), 'b': (1, 2), 'e': (0.6666666666666666, 0.8888888888888888), 'd': (1.0, 0.6666666666666666), 'f': (1.3333333333333333, 0.8888888888888888)}
        
        sage: ptsdict_2=matroids_plot_helpers.addnontripts(['a','b','c'],['d','e','f','g','h'],ptsdict)
        sage: print ptsdict_2
        {'a': (0, 0), 'c': (2, 0), 'b': (1, 2), 'e': (0.6666666666666666, 0.8888888888888888), 'd': (1.0, 0.6666666666666666), 'g': (1.0, 0.2222222222222222), 'f': (1.3333333333333333, 0.8888888888888888), 'h': (0.5555555555555555, 0.5185185185185185)}
    
    .. NOTE::

            This method does NOT do any checks.
            
    """
    tripts=[list(ptsdict[p]) for p in tripts_labels]
    pairs = [[0,1],[1,2],[0,2]]
    q = [tripts]
    num = len(nontripts_labels)
    gridpts = [[float((tripts[0][0]+tripts[1][0]+tripts[2][0])/3),float(tripts[0][1]+tripts[1][1]+tripts[2][1])/3]]
    n=0
    while n<num:
        g = trigrid(q[0])
        q.extend([[g[0],q[0][pairs[0][0]],q[0][pairs[0][1]]],[g[0],q[0][pairs[1][0]],q[0][pairs[1][1]]],[g[0],q[0][pairs[2][0]],q[0][pairs[2][1]]]])
        q.remove(q[0])
        gridpts.extend(g[1:])
        n = n + 4 
    gridpts=gridpts[0:num]
    j = 0
    for p in nontripts_labels:
        ptsdict[p]=tuple(gridpts[j])
        j = j + 1
    return ptsdict

def createline(ptsdict,ll,lineorders2=None):
    """
    Return ordered lists of co-ordinates of points to be traversed to draw a 2D line 
    INPUT:
    
    - ``ptsdict`` -- A dictionary containing keys and their (x,y) position as values.
    - ``ll`` -- A list of keys in ``ptsdict`` through which a line is to be drawn.
    - ``lineorders2``-- (optional) A list of ordered lists of keys in ``ptsdict`` such that 
    if ll is setwise same as any of these then points corresponding to values of the keys 
    will be traversed in that order thus overriding internal order deciding heuristic. 
    
    OUTPUT:
    
    A tuple containing 4 elements in this order:
    
    1. Ordered list of x-coordinates of values of keys in ``ll`` specified in ptsdict 
    2. Ordered list of y-coordinates of values of keys in ``ll`` specified in ptsdict
    3. Ordered list of interpolated x-coordinates of points through which a line can be drawn
    4. Ordered list of interpolated y-coordinates of points through which a line can be drawn    
    
    Examples::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: ptsdict={'a':(1,3),'b':(2,1),'c':(4,5),'d':(5,2)}
        sage: x,y,x_i,y_i=matroids_plot_helpers.createline(ptsdict,['a','b','c','d'])
        sage: [len(x),len(y),len(x_i),len(y_i)]
        [4, 4, 100, 100]
        sage: G = line(zip(x_i, y_i),color='black',thickness=3,zorder=1)
        sage: G+=points(zip(x, y), color='black', size=300,zorder=2)
        sage: G.show()
        
        sage: x,y,x_i,y_i=matroids_plot_helpers.createline(ptsdict,['a','b','c','d'],lineorders2=[['b','a','c','d'],['p','q','r','s']])
        sage: [len(x),len(y),len(x_i),len(y_i)]
        [4, 4, 100, 100]
        sage: G = line(zip(x_i, y_i),color='black',thickness=3,zorder=1)
        sage: G+=points(zip(x, y), color='black', size=300,zorder=2)
        sage: G.show()

        
    .. NOTE::

            This method does NOT do any checks.
            
    """
    x,lo = line_hasorder(ll,lineorders2)
    flip = False
    if x==False:
        linepts = [list(ptsdict[i]) for i in ll] # convert dictionary to list of lists
        xpts = [x[0] for x in linepts]
        ypts = [y[1] for y in linepts]
        if (float(max(xpts))-float(min(xpts))) > (float(max(ypts))-float(min(ypts))):
            sortedind = sorted(range(len(xpts)), key=lambda k: float(xpts[k]))
        else:
            sortedind = sorted(range(len(ypts)), key=lambda k: float(ypts[k]))
            flip = True
        sortedlinepts = [linepts[i] for i in sortedind]
        sortedx = [k[0] for k in sortedlinepts]
        sortedy = [k[1] for k in sortedlinepts]
    else:
        linepts = [list(ptsdict[i]) for i in lo]
        sortedx = [k[0] for k in linepts]
        sortedy = [k[1] for k in linepts]
        
    if flip==True:
        tck,u = scipy.interpolate.splprep([sortedy,sortedx],s=0.0,k=2)
        y_i,x_i = scipy.interpolate.splev(np.linspace(0,1,100),tck)
    else:
        tck,u = scipy.interpolate.splprep([sortedx,sortedy],s=0.0,k=2)
        x_i,y_i = scipy.interpolate.splev(np.linspace(0,1,100),tck)
    return sortedx,sortedy,x_i,y_i

def slp(M1):
    """
    Return simple matroid, loops and parallel elements of given matroid
    
    INPUT:
    
    - ``M1`` -- A matroid.
    
    OUTPUT:
    
    A tuple containing 3 elements in this order:
    1. Simple matroid corresponding to M1
    2. Loops of matroid ``M1``
    3. Elements that are in `M1.groundset()` but not in ground set of 1. or in 2.
    
    EXAMPLES::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: M=Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1],[0, 1, 0, 1, 0, 1, 1,0,0],[0, 0, 1, 1, 1, 0, 1,0,0]])
        sage: [M1,L,P]=matroids_plot_helpers.slp(M)
        sage: print L,P
        set([7]) set([8])
        sage: M=Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1],[0, 1, 0, 1, 0, 1, 1,0,0],[0, 0, 1, 1, 1, 0, 1,0,0]])
        sage: [M1,L,P]=matroids_plot_helpers.slp(M)
        sage: M1.is_simple()
        True
    
    .. NOTE::

            This method does NOT do any checks.
            
    """
    L = set(M1.loops())
    sg = sorted(M1.simplify().groundset())
    nP = L|set(M1.simplify().groundset())
    P = set(M1.groundset())-nP
    return [M1.delete(L|P),L,P]
    
def addlp(M,L,P,ptsdict,G=None):
    """
    Return a graphics object containing loops (in inset) and parallel elements of matroid 
    
    INPUT:
    
    - ``M`` -- A matroid.
    - ``L`` -- List of labels of elements in ``M.groundset()`` that are loops of matroid ``M``
    - ``P`` -- List of labels of elements in ``M.groundset()`` not in ``M.simplify.groundset()`` or ``L`` 
    - ``ptsdict`` -- A dictionary containing elements in ``M.groundset()`` not necessarily containing elements of ``L``
    - ``G`` -- (optional) A sage graphics object to which loops and parallel elements of matroid `M` added 
    
    OUTPUT:
    A sage graphics object containing loops and parallel elements of matroid ``M``
    
    EXAMPLES:
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: M=Matroid(ring=GF(2), matrix=[[1, 0, 0, 0, 1, 1, 1,0,1],[0, 1, 0, 1, 0, 1, 1,0,0],[0, 0, 1, 1, 1, 0, 1,0,0]])
        sage: [M1,L,P]=matroids_plot_helpers.slp(M)
        sage: G=matroids_plot_helpers.addlp(M,L,P,{0:(0,0)})
        sage: G.show(axes=False)
        sage: sage: G.show(axes=False)
    
    .. NOTE::

            This method does NOT do any checks.
            
    """
    if G == None:
        G=Graphics()
    M1 = M.simplify()
    # deal with loops
    if len(L)>0:
        loops = L
        looptext = ", ".join([str(l) for l in loops])
        rectx = -1
        recty = -1
        rectw = 0.5 + 0.4*len(loops) + 0.5 # control this based on len(loops) 
        recth = 0.6
        G += polygon2d([[rectx,recty], [rectx,recty+recth], [rectx+rectw,recty+recth], [rectx+rectw,recty]], color='black',fill=False,thickness=4)
        G += text(looptext,(rectx+0.5,recty+0.3),color='black',fontsize=13)
        G += point((rectx+0.2, recty+0.3),color='black', size=300,zorder=2)
        G += text('Loop(s)',(rectx+0.5+0.4*len(loops)+0.1,recty+0.3),fontsize=13,color='black')
    if len(P)>0:
        # create list of lists where inner lists are parallel classes 
        pcls=[]
        gnd = sorted(M1.groundset_list())
        for g in gnd:
            pcl = [g]
            for p in P:
                if M.rank([g,p])==1:
                    pcl.extend([p])
            pcls.append(pcl)
        for pcl in pcls:
            if len(pcl)>1:
                basept = list(ptsdict[pcl[0]])
                if len(pcl)<=2:
                    # add side by side
                    ptsdict[pcl[1]] = (basept[0],basept[1]-0.13)
                    G += points(zip([basept[0]], [basept[1]-0.13]),color='black', size=300,zorder=2)
                    G += text(pcl[1],(float(basept[0]), float(basept[1])-0.13),color='white',fontsize=13)
                else:
                    # add in a bracket
                    G += text('{ '+", ".join(sorted([str(kk) for kk in pcl[1:]]))+' }',(float(basept[0]), float(basept[1]-0.2)-0.034),color='black',fontsize=13)
    return G
      
def line_hasorder(l,lodrs=None):
    """
    Determine if and order is specified for a line
    
    INPUT:
    
    - ``l`` -- A line containing specified as a list of labels
    - ``lordrs`` -- (optional) A list of lists each specifying an order on the subset of labels corresponding to that list
     
    OUTPUT:
    
    A tuple containing 3 elements in this order:
    1. A boolean indicating whether there is any list in ``lordrs`` that is setwise equal to ``l``
    2. A list specifying an order on ``set(l)`` if 1. is True, otherwise an empty list
    
    EXAMPLES::
        
        
        sage: matroids_plot_helpers.line_hasorder(['a','b','c','d'],[['a','c','d','b'],['p','q','r']])
        (True, ['a', 'c', 'd', 'b'])
        sage: matroids_plot_helpers.line_hasorder(['a','b','c','d'],[['p','q','r'],['l','m','n','o']])
        (False, [])
        
    .. NOTE::

            This method does NOT do any checks.
            
    """
    if lodrs!=None:
        if len(lodrs) > 0:
            for i in lodrs:
                if Set(i)==Set(l):
                    return True,i
    return False,[]    
    
def geomrep(M1,B1=None,lineorders1=None):
    """
    Return a sage graphics object containing geometric representation of matroid M1
    
    INPUT:
    
    - ``M1`` -- A matroid.
    - ``B1`` -- (optional) A list of labels in ``M1.groundset()`` that correspond to a basis 
    of ``M1`` and will be placed as vertices of the triangle in the geometric representation of ``M1``
    - ``lineorders1`` -- (optional) A list of ordered lists of elements of ``M1.grondset()`` such that 
    if a line in geometric representation is setwise same as any of these then points contained will be
    traversed in that order thus overriding internal order deciding heuristic.
    
    OUTPUT:
    
    A sage graphics object of type <class 'sage.plot.graphics.Graphics'> that 
    corresponds to the geometric representation of the matroid
    
    EXAMPLES::
        
        sage: from sage.matroids import matroids_plot_helpers
        sage: M=matroids.named_matroids.P7()
        sage: G=matroids_plot_helpers.geomrep(M)
        sage: G.show(xmin=-2, xmax=3, ymin=-2, ymax=3)
        
        sage: M=matroids.named_matroids.P7()
        sage: G=matroids_plot_helpers.geomrep(M,lineorders1=[['f','e','d']])
        sage: G.show(xmin=-2, xmax=3, ymin=-2, ymax=3)

    
    """
    if B1 == None:
        B1 = list(M1.basis())
    G = Graphics()
    # create lists of loops and parallel elements and simplify given matroid
    [M,L,P] = slp(M1)
    if M.rank()==1:
        pts = {}
        gnd = sorted(M.groundset())
        pts[gnd[0]]=(1,float(2)/3)
        G += point((1, float(2)/3),size=300,zorder=2)
        pts2 = pts
    elif M.rank() == 2:
        pts2 = {}
        pts2[B1[0]] = (0,0)
        pts2[B1[1]] = (2,0)
        nB1=list(set(M.groundset_list())-set(B1))
        bline = []
        for j in nB1:
            if M.is_dependent([j,B1[0],B1[1]]):
                bline.append(j)
        interval = len(bline)+1
        lpt = list(pts2[B1[0]])
        rpt = list(pts2[B1[1]])
        for k in range(len(bline)):
            cc = (float(1)/interval)*(k+1)
            pts2[bline[k]] = (cc*lpt[0]+(1-cc)*rpt[0],cc*lpt[1]+(1-cc)*rpt[1])
        bline.extend(B1)
        ptsx,ptsy,x_i,y_i = createline(pts2,bline,lineorders1)
        G += line(zip(x_i, y_i),color='black',thickness=3,zorder=1)
        allpts = [list(pts2[i]) for i in M.groundset()]
        xpts = [float(k[0]) for  k in allpts]
        ypts = [float(k[1]) for  k in allpts]
        G += points(zip(xpts, ypts), color='black', size=300,zorder=2)
        for i in pts2:
            pt = list(pts2[i])
            G += text(i,(float(pt[0]), float(pt[1])), color='white',fontsize=13)    
    else:
        pts,trilines,nontripts,curvedlines = initial_triangle(M1,B1,list(set(M.groundset())-set(B1)), list(set(L)|set(P))) #[i for i in sorted(M.groundset()) if i not in B1])
        pts2= addnontripts([B1[0],B1[1],B1[2]],nontripts,pts)
        trilines.extend(curvedlines)
        j = 0
        for ll in trilines:
            if len(ll)>=3:
                ptsx,ptsy,x_i,y_i=createline(pts2,ll,lineorders1)
                G += line(zip(x_i,y_i),color='black',thickness=3,zorder=1)
        allpts = [list(pts2[i]) for i in M.groundset()]
        xpts = [float(k[0]) for  k in allpts]
        ypts = [float(k[1]) for  k in allpts]
        G += points(zip(xpts,ypts),color='black', size=300,zorder=2)
        for i in pts2:
            pt = list(pts2[i])
            G += text(i,(float(pt[0]), float(pt[1])),color='white',fontsize=13)
    #deal with loops and parallel elements 
    G = addlp(M1,L,P,pts2,G)
    G.axes(False)
    return G
