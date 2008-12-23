r"""
Covering designs: coverings of $t$-element subsets of a $v$-set by $k$-sets

A $(v,k,t)$ covering design $C$ is an incidence structure consisting of a
set of points $P$ of order $v$, and a set of blocks $B$, where each
block contains $k$ points of $P$.  Every $t$-element subset of $P$
must be contained in at least one block.

If every $t$-set is contained in exactly one block of $C$, then we
have a block design.  Following the block design implementation, the
standard representation of a covering design uses $P = [0,1,..., v-1]$.

There is an online database of the best known covering designs, the La
Jolla Covering Repository (LJCR), at [1].  This is maintained with an SQL
database, also in Sage, but since it changes frequently, and is over
60MB, that code is not included here.  A user may get individual
coverings from the LJCR using best_known_covering_design_www.

In addition to the parameters and incidence structure for a covering
design from this database, we include extra information:

\begin{itemize}
\item Best known lower bound on the size of a (v,k,t)-covering design
\item Name of the person(s) who produced the design
\item Method of construction used
\item Date when the design was added to the database
\end{itemize}

REFERENCES
  [1] La Jolla Covering Repository, http://www.ccrwest.org/cover.html
  [2] Coverings, by Daniel Gordon and Douglas Stinson,
      http://www.ccrwest.org/gordon/hcd.pdf, in Handbook of Combinatorial Designs


AUTHORS:
    -- Daniel M. Gordon (2008-12-22): initial version

"""

#*****************************************************************************
#      Copyright (C) 2008 Daniel M. Gordon <dmgordo@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

import types
import urllib
from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial
from sage.combinat.combination import Combinations
from sage.rings.real_double import RDF
from sage.combinat.designs.incidence_structures import IncidenceStructure

###################### covering design functions ##############################


def schonheim(v,k,t):
    r"""
    Schonheim lower bound for size of covering design


    INPUT:
        v -- integer, size of point set
        k -- integer, cardinality of each block
        t -- integer, cardinality of sets being covered

    OUTPUT:
        The Schonheim lower bound for such a covering design's size:
        $C(v,k,t) \leq \lceil(\frac{v}{k}  \lceil \frac{v-1}{k-1} \cdots \lceil \frac{v-t+1}{k-t+1} \rceil \cdots \rceil \rceil$

    EXAMPLES:
        sage: from sage.combinat.designs.covering_design import schonheim
        sage: schonheim(10,3,2)
        17
        sage: schonheim(32,16,8)
        930
    """
    bound = 1.0
    for i in range(t-1,-1,-1):
        bound = (bound*RDF(v-i)/RDF(k-i)).ceiling()

    return bound


def trivial_covering_design(v,k,t):
    r"""
    Construct a trivial covering design.

    INPUT:
        v -- integer, size of point set
        k -- integer, cardinality of each block
        t -- integer, cardinality of sets being covered

    OUTPUT:
        (v,k,t) covering design

    EXAMPLE:
        sage: C = trivial_covering_design(8,3,1)
        sage: C.show()
        C(8,3,1) = 3
        Method: Trivial
        0   1   2
        0   6   7
        3   4   5
        sage: C = trivial_covering_design(5,3,2)
        sage: C.show()
        4 <= C(5,3,2) <= 10
        Method: Trivial
        0   1   2
        0   1   3
        0   1   4
        0   2   3
        0   2   4
        0   3   4
        1   2   3
        1   2   4
        1   3   4
        2   3   4

    NOTES:
        Cases are:
        \begin{enumerate}
        \item $t=0$: This could be empty, but it's a useful convention to
        have one block (which is empty if $k=0$).
        \item $t=1$: This contains $\lceil v/k \rceil$ blocks:
        $[0,...,k-1]$,[k,...,2k-1],...$.  The last block
        wraps around if $k$ does not divide $v$.
        \item anything else:  Just use every $k$-subset of $[0,1,...,v-1]$.
        \end{enumerate}

    """
    if t==0:     #single block [0,...,k-1]
        blk=[]
        for i in range(k):
            blk.append(i)
        return CoveringDesign(v,k,t,1,range(v),[blk],1,"Trivial")

    if t==1:     #blocks [0,...,k-1],[k,...,2k-1],...
        size = (RDF(v)/RDF(k)).ceiling()
        blocks=[]
        for i in range(size-1):
            blk=[]
            for j in range(i*k,(i+1)*k):
                blk.append(j)
            blocks.append(blk)

        blk=[]   # last block: if k does not divide v, wrap around
        for j in range((size-1)*k,v):
            blk.append(j)
        for j in range(k-len(blk)):
            blk.append(j)
        blk.sort()
        blocks.append(blk)
        return CoveringDesign(v,k,t,size,range(v),blocks,size,"Trivial")

    # default case, all k-subsets
    return CoveringDesign(v,k,t,binomial(v,k),range(v),Combinations(range(v),k),schonheim(v,k,t),"Trivial")


class CoveringDesign(SageObject):
    """
    Covering design.

    INPUT:
      data -- parameters (v,k,t)
              size
              incidence structure (points and blocks) -- default points are $[0,...v-1]$
              lower bound for such a design
              database information: creator, method, timestamp
    """

    def __init__(self, v=0, k=0, t=0, size=0, points=[], blocks=[], low_bd=0, method='', created_by ='',timestamp=''):
        """
        EXAMPLES:
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.show()
            C(7,3,2) = 7
            Method: Projective Plane
            0   1   2
            0   3   4
            0   5   6
            1   3   5
            1   4   6
            2   3   6
            2   4   5
        """
        self.v = v
        self.k = k
        self.t = t
        self.size = size
        if low_bd > 0:
            self.low_bd = low_bd
        else:
            self.low_bd = schonheim(v,k,t)
        self.method = method
        self.created_by = created_by
        self.timestamp = timestamp
        self.incidence_structure = IncidenceStructure(points, blocks)


    def show(self, max_size=100):
        """
        Displays a covering design.

        INPUT:
            max_size -- maximum number of blocks (to avoid trying to show huge ones)

        OUTPUT:
            covering design parameters and blocks, in a readable form

        EXAMPLES:
            sage: C=CoveringDesign(5,4,3,4,range(5),[[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4]],4, 'Lexicographic Covering')
            sage: C.show()
            C(5,4,3) = 4
            Method: Lexicographic Covering
            0   1   2   3
            0   1   2   4
            0   1   3   4
            0   2   3   4
        """
        if self.size == self.low_bd:    # check if the covering is known to be optimal
            print 'C(%d,%d,%d) = %d'%(self.v,self.k,self.t,self.size)
        else:
            print '%d <= C(%d,%d,%d) <= %d'%(self.low_bd,self.v,self.k,self.t,self.size);
        if self.created_by != '':
            print 'Created by: %s'%(self.created_by)
        if self.method != '':
            print 'Method: %s'%(self.method)
        if self.timestamp != '':
            print 'Submitted on: %s'%(self.timestamp)
        for i in range(min(max_size, self.size)):
            for j in range(self.k):
                print self.incidence_structure.blocks()[i][j]," ",
            print "\n",
        if(max_size < self.size):   # if not showing all blocks, indicate it with ellipsis
            print "..."



    def is_covering(self):
        """
        Checks that all t-sets are in fact covered.

        INPUT:
            putative covering design

        OUTPUT:
            True if all t-sets are in at least one block


        NOTES:
            This is very slow and wasteful of memory.  A faster cython
            version will be added soon.

        EXAMPLES:
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.is_covering()
            True
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 6]],0, 'not a covering')   # last block altered
            sage: C.is_covering()
            False
        """
        v = self.v
        k = self.k
        t = self.t
        Svt = Combinations(range(v),t)
        Skt = Combinations(range(k),t)
        tset = {}       # tables of t-sets: False = uncovered, True = covered
        for i in Svt:
            tset[tuple(i)] = False

        # mark all t-sets covered by each block
        for a in self.incidence_structure.blocks():
            for z in Skt:
                y = [a[x] for x in z]
                tset[tuple(y)] = True

        for i in Svt:
            if tset[tuple(i)] == False:     # uncovered
                return False

        return True                  # everything was covered






def best_known_covering_design_www(v, k, t, verbose=False):
    r"""
    Gives the best known (v,k,t) covering design, using database at
         \url{http://www.ccrwest.org/}.

    INPUT:
        v -- integer, the size of the point set for the design
        k -- integer, the number of points per block
        t -- integer, the size of sets covered by the blocks
        verbose -- bool (default=False), print verbose message

    OUTPUT:
        CoveringDesign -- (v,k,t) covering design with smallest number of blocks

    EXAMPLES:
        sage: C = best_known_covering_design_www(7, 3, 2)   # requires internet, optional
        sage: C.show()                                      # requires internet, optional
        C(7,3,2) = 7
        Method: lex covering
        Submitted on: 1996-12-01 00:00:00
        0   1   2
        0   3   4
        0   5   6
        1   3   5
        1   4   6
        2   3   6
        2   4   5

    This function raises a ValueError if the (v,k,t) parameters are
    not found in the database.
    """
    v = int(v)
    k = int(k)
    t = int(t)

    param = ("?v=%s&k=%s&t=%s"%(v,k,t))

    url = "http://www.ccrwest.org/cover/get_cover.php"+param
    if verbose:
        print "Looking up the bounds at %s"%url
    f = urllib.urlopen(url)
    s = f.read()
    f.close()

    if 'covering not in database' in s:   #not found
        str = "no (%d,%d,%d) covering design in database\n"%(v,k,t)
        raise ValueError, str

    return sage_eval(s)

