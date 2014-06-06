import scipy
import scipy.interpolate

def initial_triangle(M,B1,nB1,lps):
    """takes a matroid of rank=3, a basis B1, non-basis elements that are not loops or parallel elements nB1 and loops lps.
    returns a dictionary of points with M.rank() elements from B1 and those spanned by any two-subset of basis as keys and 
    their cartesian co-ordinates in R^2 as values along with non-triangle points and remaining lines that will be (mostly) 
    drawn curved
    """ 
    tripts=[(0, 0),(1, 2),(2, 0)]
    pts={}
    j=0
    for i in B1:
        pts[i]=tripts[j]
        j=j+1
    #nB1=[j for j in gnd if j not in B1]
    pairs=[[0, 1],[1,2],[0,2]]
    L1=[]
    L2=[]
    L3=[]
    for i in nB1:
        if M.is_dependent([i,B1[pairs[0][0]],B1[pairs[0][1]]]):
            L1.append(i)
        elif M.is_dependent([i,B1[pairs[1][0]],B1[pairs[1][1]]]):
            # Add to L2  
            L2.append(i)
        elif M.is_dependent([i,B1[pairs[2][0]],B1[pairs[2][1]]]):
            # Add to L3 
            L3.append(i)
    L=[L1,L2,L3]
    lines=[] #the list of lines
    for i in range(1,len(L)+1):
        lines.append([B1[pairs[i-1][0]]])
        lines[i-1].extend(L[i-1])
        lines[i-1].extend([B1[pairs[i-1][1]]])
    # place triangle and L1,L2,L3
    for i in L:# loop over megalist
        interval=1/(len(i)+1)
        pt1=list(tripts[pairs[L.index(i)][0]])
        pt2=list(tripts[pairs[L.index(i)][1]])
        for j in range(1,len(i)+1): 
            #Loop over L1,L2,L3
            cc= interval*j
            pts[i[j-1]]=(cc*pt1[0]+(1-cc)*pt2[0],cc*pt1[1]+(1-cc)*pt2[1])
    trilines=[list(set(x)) for x in lines if len(x)>=3]
    curvedlines = [list(set(list(x)).difference(set(lps))) for x in M.flats(2) if set(list(x)) not in trilines if len(list(x))>=3]
    nontripts= [i for i in nB1 if i not in pts.keys()] 
    return pts,trilines,nontripts,curvedlines

def trigrid(tripts):
    """tripts is a list of 3 points(each a 2-tuple) in R^2. returns 4 points that are in the convex hull of tripts 
        and are arranged in a grid inside the triangle 
    """ 
    n=0
    pairs=[[0,1],[1,2],[0,2]]
    cpt=list(((tripts[0][0]+tripts[1][0]+tripts[2][0])/3,(tripts[0][1]+tripts[1][1]+tripts[2][1])/3))
    grid=[cpt]
    for p in pairs:
        pt=list(((tripts[p[0]][0]+tripts[p[1]][0]+cpt[0])/3,(tripts[p[0]][1]+tripts[p[1]][1]+cpt[1])/3))
        grid.append(pt)
    return grid
    
def addnontripts(M,ptsdict,tripts,nontripts):
    pairs=[[0,1],[1,2],[0,2]]
    q=[tripts]
    num=len(nontripts)
    gridpts=[[(tripts[0][0]+tripts[1][0]+tripts[2][0])/3,(tripts[0][1]+tripts[1][1]+tripts[2][1])/3]]
    n=0
    while n<num:
        g=trigrid(q[0])
        q.extend([[g[0],q[0][pairs[0][0]],q[0][pairs[0][1]]],[g[0],q[0][pairs[1][0]],q[0][pairs[1][1]]],[g[0],q[0][pairs[2][0]],q[0][pairs[2][1]]]])
        q.remove(q[0])
        gridpts.extend(g[1:])
        n=n+4 
    gridpts=gridpts[0:num]
    j=0
    for p in nontripts:
        ptsdict[p]=tuple(gridpts[j])
        j=j+1
    return ptsdict

def createline(ll,ptsdict,lineorders=None):
    """
    Flips x-y if necessary and returns x,y, x_i, y_i lists to be passed to plot function 
    """
    x,lo= line_hasorder(ll,lineorders)
    flip = False
    if x==False:
        linepts=[list(ptsdict[i]) for i in ll] # convert dictionary to list of lists
        xpts=[x[0] for x in linepts]
        ypts=[y[1] for y in linepts]
        if (float(max(xpts))-float(min(xpts))) > (float(max(ypts))-float(min(ypts))):
            sortedind= sorted(range(len(xpts)), key=lambda k: float(xpts[k]))
        else:
            sortedind= sorted(range(len(ypts)), key=lambda k: float(ypts[k]))
            flip = True
        sortedlinepts=[linepts[i] for i in sortedind]
        sortedx=[k[0] for k in sortedlinepts]
        sortedy=[k[1] for k in sortedlinepts]
    else:
        linepts=[list(ptsdict[i]) for i in lo]
        sortedx=[k[0] for k in linepts]
        sortedy=[k[1] for k in linepts]
        
    if flip==True:
        tck,u=scipy.interpolate.splprep([sortedy,sortedx],s=0.0,k=2)
        y_i,x_i= scipy.interpolate.splev(np.linspace(0,1,100),tck)
    else:
        tck,u=scipy.interpolate.splprep([sortedx,sortedy],s=0.0,k=2)
        x_i,y_i= scipy.interpolate.splev(np.linspace(0,1,100),tck)
    return sortedx,sortedy,x_i,y_i

def slp(M1):
    """Takes a matroid M1 and returns it corresponding simple matroid, loops and parallel elements
    """
    L=set(M1.loops())
    sg=sorted(M1.simplify().groundset())
    nP=L|set(M1.simplify().groundset())
    P= set(M1.groundset())-nP
    return [M1.delete(L|P),L,P]
    
def addlp(G,M,ptsdict,L,P):
    """Takes a graphics object G, points dictionary, loops and parallel elements and adds loops and 
    parallel elements to graphics object and returns the modified graphics object
    """
    M1=M.simplify()
    # deal with loops
    if len(L)>0:
        loops=L
        looptext=", ".join([str(l) for l in loops])
        rectx=-1
        recty=-1
        rectw= 0.5+0.4*len(loops)+ 0.5 # control this based on len(loops) 
        recth=0.6
        G+=polygon2d([[rectx,recty], [rectx,recty+recth], [rectx+rectw,recty+recth], [rectx+rectw,recty]], color='black',fill=False,thickness=4)
        G+=text(looptext,(rectx+0.5,recty+0.3),color='black',fontsize=13)
        G+=point((rectx+0.2, recty+0.3),color='black', size=300,zorder=2)
        G+=text('Loop(s)',(rectx+0.5+0.4*len(loops)+0.1,recty+0.3),fontsize=13,color='black')
    if len(P)>0:
        # create list of lists where inner lists are parallel classes 
        pcls=[]
        gnd=sorted(M1.groundset_list())
        for g in gnd:
            pcl=[g]
            for p in P:
                if M.rank([g,p])==1:
                    pcl.extend([p])
            pcls.append(pcl)
        # add parallel elements to axis object    
        for pcl in pcls:
            if len(pcl)>1:
                basept=list(ptsdict[pcl[0]])
                if len(pcl)<=2:
                    # add side by side
                    ptsdict[pcl[1]]=(basept[0],basept[1]-0.13)
                    G+=points(zip([basept[0]], [basept[1]-0.13]),color='black', size=300,zorder=2)
                    G+=text(pcl[1],(float(basept[0]), float(basept[1])-0.13),color='white',fontsize=13)
                else:
                    # add in a bracket
                    G+=text('{ '+", ".join(sorted([str(kk) for kk in pcl[1:]]))+' }',(float(basept[0]), float(basept[1]-0.2)-0.034),color='black',fontsize=13)
    return G
      
def line_hasorder(l,lodrs=None):
    if lodrs!=None:
        if len(lodrs) > 0:
            for i in lodrs:
                if Set(i)==Set(l):
                                        return True,i
    return False,[]    
