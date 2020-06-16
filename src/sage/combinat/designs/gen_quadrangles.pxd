def generalised_quadrangle_with_spread(const int s, const int t, existence=False, check=True):
    r"""
    Returns a pair (GQ,S) s.t. GQ is a generalised quadrangle of order (s,t) and S is a spread of GQ
    """
    if s < 1 or t < 1:
        if existence: return False
        raise RuntimeError("No GQ of order ({},{}) exists".format(s,t))

    if s == 1 and t == 1:#we have a square
        if existence: return True
        D = IncidenceStructure([[0,1],[1,2],[2,3],[3,0]])
        return (D,[[0,1],[2,3]])

    if is_prime_power(s) and t == s*s:
        if existence: return True
        (GQ,S) = dual_GQ_ovoid(*generalised_quadrangle_hermitian(s))
        if check:
            if not is_GQ_with_spread(GQ,S,s,t):
                raise RuntimeError("Sage built a wrong GQ with spread")
        return (GQ,S)

    if existence: return Unknown
    raise RuntimeError("Sage can't build a GQ of order ({},{}) with a spread".format(s,t))

def is_GQ_with_spread(GQ,S,const int s, const int t):
    r"""
    Checks if GQ is a generalised quadrangle of order (s,t) and
    checks that S is a spred of GQ
    """
    res = GQ.is_generalised_quadrangle(parameters=True)
    if res is False or res[0] != s or res[1] != t:
        return False

    #check spread
    points = set(GQ.ground_set())
    for line in S:
        if not points.issuperset(line):
            return False
        points = points.difference(line)
        
    if points:
        return False
    
    return True
        
def dual_GQ_ovoid(GQ,O):
    r"""
    Computes the dual of GQ and returns the image of O under the dual map
    """
    #we compute the dual of GQ and of O

    #GQ.ground_set()[i] becomes newBlocks[i]
    #GQ.blocks()[i] becomes i
    newBlocks = [ [] for _ in range(GQ.num_points())]
    pointsToInt = { p: i for i,p in enumerate(GQ.ground_set()) }

    for i,b in enumerate(GQ.blocks()):
        for p in b:
            newBlocks[pointsToInt[p]].append(i)

    S = [ newBlocks[pointsToInt[p]] for p in O]

    D = IncidenceStructure(newBlocks)
    return (D,S)
    
def generalised_quadrangle_hermitian(const int q):
    r"""
    Construct the generalised quadrangle H(3,q^2) with an ovoid
    The GQ has order (q^2,q)
    """

    GU = libgap.GU(4,q)
    H = libgap.InvariantSesquilinearForm(GU)["matrix"]
    Fq = libgap.GF(q*q)
    zero = libgap.Zero(Fq)
    one = libgap.One(Fq)
    V = libgap.FullRowSpace(Fq,4)

    e1 = [one,zero,zero,zero] #isotropic point
    assert( e1*H*e1 == zero, "e1 not isotropic")

    points = list(libgap.Orbit(GU,e1,libgap.OnLines)) #all isotropic points
    pointInt = { x:(i+1) for i,x in enumerate(points) } #+1 because GAP starts at 1
    #points is the hermitian variety

    GUp = libgap.Action(GU, points, libgap.OnLines)#GU as permutation group of points

    e2 = [zero,one,zero,zero]
    #we have totally isotropic line
    line = V.Subspace([e1,e2])
    lineAsPoints = [libgap.Elements(libgap.Basis(b))[0] for b in libgap.Elements(line.Subspaces(1)) ]
    line = libgap.Set([ pointInt[p] for p in lineAsPoints ])

    lines = libgap.Orbit(GUp, line, libgap.OnSets)#all isotropic lines

    #to find ovoid, we embed H(3,q^2) in H(4,q^2)
    #then embedding is (a,b,c,d) -> (a,b,0,c,d) [so we preserve isotropicity]
    #then find a point in the latter and not in the former
    #this point will be collinear in H(3,q^2) to all (and only) the points in a ovoid
    W = libgap.FullRowSpace(Fq,5)
    J = [ [0,0,0,0,1],[0,0,0,1,0],[0,0,1,0,0],[0,1,0,0,0],[1,0,0,0,0]]
    J = libgap(J)
    if q%2 == 1:
        (p,k) = is_prime_power(q,get_data=True)
        a = (p-1)// 2
        aGap = zero
        for i in range(a): aGap += one
        p = [zero,one,one,aGap,zero]
    else:
        a = libgap.PrimitiveRoot(Fq)**(q-1)
        p = [zero,one,a+one,a,zero]
        
    #now p is a point of H(4,q^2)

    #p' is collinear to p iff p'Jp^q = 0
    #note that p'Jp^q = bx^q + c where p' =(a,b,0,c,d) and p=(0,1,1,x,0)
    #hece we have points (0,0,0,1); (0,1,c,a) for any a iff c^q+c = 0 (c = -x^q)
    #and (1,b,c,x) for any x and any b (c= -bx^q) iff it is an iso point
    #so we need only q^2 (for last case) +1 (for 2nd case) checks 
    ovoid = []
    xq = p[3]**q
    for p2 in points:
        if p2[1]*xq+p2[2] == zero: #p collinear to newP2
            ovoid.append(libgap(pointInt[p2]))

    D = IncidenceStructure(lines)
    return (D,ovoid)
