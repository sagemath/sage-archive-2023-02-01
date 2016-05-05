r"""
Covering designs: coverings of `t`-element subsets of a `v`-set by `k`-sets

A `(v,k,t)` covering design `C` is an incidence structure consisting of a
set of points `P` of order `v`, and a set of blocks `B`, where each
block contains `k` points of `P`.  Every `t`-element subset of `P`
must be contained in at least one block.

If every `t`-set is contained in exactly one block of `C`, then we
have a block design.  Following the block design implementation, the
standard representation of a covering design uses `P = [0,1,..., v-1]`.

In addition to the parameters and incidence structure for a covering
design from this database, we include extra information:

* Best known lower bound on the size of a `(v,k,t)`-covering design
* Name of the person(s) who produced the design
* Method of construction used
* Date when the design was added to the database

REFERENCES:

.. [1] La Jolla Covering Repository,
  http://www.ccrwest.org/cover.html

.. [2] Coverings,
  Daniel Gordon and Douglas Stinson,
  http://www.ccrwest.org/gordon/hcd.pdf
  from the Handbook of Combinatorial Designs

AUTHORS:

    -- Daniel M. Gordon (2008-12-22): initial version

Classes and methods
-------------------
"""

#*****************************************************************************
#      Copyright (C) 2008 Daniel M. Gordon <dmgordo@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject
from sage.rings.rational import Rational
from sage.arith.all import binomial
from sage.combinat.combination import Combinations
from sage.combinat.designs.incidence_structures import IncidenceStructure

###################### covering design functions ##############################


def schonheim(v,k,t):
    r"""
    Schonheim lower bound for size of covering design

    INPUT:

    - ``v`` -- integer, size of point set
    - ``k`` -- integer, cardinality of each block
    - ``t`` -- integer, cardinality of sets being covered

    OUTPUT:

    The Schonheim lower bound for such a covering design's size:
    `C(v,k,t) \leq \lceil(\frac{v}{k}  \lceil \frac{v-1}{k-1} \cdots \lceil \frac{v-t+1}{k-t+1} \rceil \cdots \rceil \rceil`

    EXAMPLES::

        sage: from sage.combinat.designs.covering_design import schonheim
        sage: schonheim(10,3,2)
        17
        sage: schonheim(32,16,8)
        930
    """
    bound = 1
    for i in range(t-1,-1,-1):
        bound = Rational((bound*(v-i),k-i)).ceil()

    return bound


def trivial_covering_design(v,k,t):
    r"""
    Construct a trivial covering design.

    INPUT:

    - ``v`` -- integer, size of point set
    - ``k`` -- integer, cardinality of each block
    - ``t`` -- integer, cardinality of sets being covered

    OUTPUT:

    `(v,k,t)` covering design

    EXAMPLE::

        sage: C = trivial_covering_design(8,3,1)
        sage: print C
        C(8,3,1) = 3
        Method: Trivial
        0   1   2
        0   6   7
        3   4   5
        sage: C = trivial_covering_design(5,3,2)
        sage: print C
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

    * `t=0`: This could be empty, but it's a useful convention to have
      one block (which is empty if $k=0$).

    * `t=1` : This contains `\lceil v/k \rceil` blocks:
      `[0,...,k-1],[k,...,2k-1],...`.  The last block wraps around if
      `k` does not divide `v`.

    * anything else: Just use every `k`-subset of `[0,1,...,v-1]`.

    """
    if t==0:     #single block [0,...,k-1]
        blk=[]
        for i in range(k):
            blk.append(i)
        return CoveringDesign(v,k,t,1,range(v),[blk],1,"Trivial")

    if t==1:     #blocks [0,...,k-1],[k,...,2k-1],...
        size = Rational((v,k)).ceil()
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

    - ``v,k,t`` -- integer parameters of the covering design.

    - ``size`` (integer)

    - ``points`` -- list of points (default points are `[0,...v-1]`).

    - ``blocks``

    - ``low_bd`` (integer) -- lower bound for such a design

    - ``method, creator, timestamp`` -- database information.
    """

    def __init__(self, v=0, k=0, t=0, size=0, points=[], blocks=[], low_bd=0, method='', creator ='',timestamp=''):
        """
        EXAMPLES::

            sage: C=CoveringDesign(5,4,3,4,range(5),[[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4]],4, 'Lexicographic Covering')
            sage: print C
            C(5,4,3) = 4
            Method: Lexicographic Covering
            0   1   2   3
            0   1   2   4
            0   1   3   4
            0   2   3   4
        """
        self.__v = v
        self.__k = k
        self.__t = t
        self.__size = size
        if low_bd > 0:
            self.__low_bd = low_bd
        else:
            self.__low_bd = schonheim(v,k,t)
        self.__method = method
        self.__creator = creator
        self.__timestamp = timestamp
        self.__incidence_structure = IncidenceStructure(points, blocks)


    def __repr__(self):
        """
        A print method, giving the parameters and any other
        information about the covering (but not the blocks).

        EXAMPLES::

            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C
            (7,3,2)-covering design of size 7
            Lower bound: 7
            Method: Projective Plane
        """
        repr = '(%d,%d,%d)-covering design of size %d\n' % (self.__v,
                                                            self.__k,
                                                            self.__t,
                                                            self.__size)
        repr += 'Lower bound: %d\n' % (self.__low_bd)
        if self.__creator != '':
            repr += 'Created by: %s\n' % (self.__creator)
        if self.__method != '':
            repr += 'Method: %s\n' % (self.__method)
        if self.__timestamp != '':
            repr += 'Submitted on: %s\n' % (self.__timestamp)

        return repr

    def __str__(self):
        """
        A print method, displaying a covering design's parameters and blocks.

        OUTPUT:

        covering design parameters and blocks, in a readable form

        EXAMPLES::

            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: print C
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
        if self.__size == self.__low_bd:    # check if the covering is known to be optimal
            repr =  'C(%d,%d,%d) = %d\n'%(self.__v,self.__k,self.__t,self.__size)
        else:
            repr = '%d <= C(%d,%d,%d) <= %d\n'%(self.__low_bd,self.__v,self.__k,self.__t,self.__size);
        if self.__creator != '':
            repr += 'Created by: %s\n'%(self.__creator)
        if self.__method != '':
            repr += 'Method: %s\n'%(self.__method)
        if self.__timestamp != '':
            repr += 'Submitted on: %s\n'%(self.__timestamp)
        for i in range(self.__size):
            for j in range(self.__k):
                repr = repr + str(self.__incidence_structure.blocks()[i][j]) + '  '
            repr += '\n'

        return repr

    def is_covering(self):
        """
        Checks that all `t`-sets are in fact covered by the blocks of
        ``self``

        .. NOTE::

            This is very slow and wasteful of memory.

        EXAMPLES::

            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.is_covering()
            True
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 6]],0, 'not a covering')   # last block altered
            sage: C.is_covering()
            False
        """
        v = self.__v
        k = self.__k
        t = self.__t
        Svt = Combinations(range(v),t)
        Skt = Combinations(range(k),t)
        tset = {}       # tables of t-sets: False = uncovered, True = covered
        for i in Svt:
            tset[tuple(i)] = False

        # mark all t-sets covered by each block
        for a in self.__incidence_structure.blocks():
            for z in Skt:
                y = [a[x] for x in z]
                tset[tuple(y)] = True

        for i in Svt:
            if not tset[tuple(i)]:     # uncovered
                return False

        return True                  # everything was covered

    def v(self):
        """
        Return `v`, the number of points in the covering design.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.v()
            7
        """
        return self.__v


    def k(self):
        """
        Return `k`, the size of blocks of the covering design

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.k()
            3
        """
        return self.__k


    def t(self):
        """
        Return `t`, the size of sets which must be covered by the
        blocks of the covering design

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.t()
            2
        """
        return self.__t

    def size(self):
        """
        Return the number of blocks in the covering design

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.size()
            7
        """
        return self.__size


    def low_bd(self):
        """
        Return a lower bound for the number of blocks a covering
        design with these parameters could have.

        Typically this is the Schonheim bound, but for some parameters
        better bounds have been shown.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.low_bd()
            7
        """
        return self.__low_bd

    def method(self):
        """
        Return the method used to create the covering design
        This field is optional, and is used in a database to give information about how coverings were constructed

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: C.method()
            'Projective Plane'
        """
        return self.__method

    def creator(self):
        """
        Return the creator of the covering design

        This field is optional, and is used in a database to give
        attribution for the covering design It can refer to the person
        who submitted it, or who originally gave a construction

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane','Gino Fano')
            sage: C.creator()
            'Gino Fano'
        """
        return self.__creator

    def timestamp(self):
        """
        Return the time that the covering was submitted to the database

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane','Gino Fano','1892-01-01 00:00:00')
            sage: C.timestamp()  #Fano had an article in 1892, I don't know the date it appeared
            '1892-01-01 00:00:00'
        """
        return self.__timestamp


    def incidence_structure(self):
        """
        Return the incidence structure of a covering design, without all the extra parameters.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_design import CoveringDesign
            sage: C=CoveringDesign(7,3,2,7,range(7),[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]],0, 'Projective Plane')
            sage: D = C.incidence_structure()
            sage: D.ground_set()
            [0, 1, 2, 3, 4, 5, 6]
            sage: D.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]

        """
        return self.__incidence_structure

def best_known_covering_design_www(v, k, t, verbose=False):
    r"""
    Gives the best known `(v,k,t)` covering design, using the database
    available at `<http://www.ccrwest.org/>`_

    INPUeT:

    - ``v`` -- integer, the size of the point set for the design
    - ``k`` -- integer, the number of points per block
    - ``t`` -- integer, the size of sets covered by the blocks
    - ``verbose`` -- bool (default=``False``), print verbose message

    OUTPUT:

    A :class:`CoveringDesign` object representing the ``(v,k,t)``-covering
    design with smallest number of blocks available in the database.

    EXAMPLES::

        sage: from sage.combinat.designs.covering_design import best_known_covering_design_www
        sage: C = best_known_covering_design_www(7, 3, 2)   # optional - internet
        sage: print C                                       # optional - internet
        C(7,3,2) = 7
        Method: lex covering
        Submitted on: 1996-12-01 00:00:00
        0  1  2
        0  3  4
        0  5  6
        1  3  5
        1  4  6
        2  3  6
        2  4  5

    This function raises a ValueError if the ``(v,k,t)`` parameters are not
    found in the database.
    """
    # import compatible with py2 and py3
    from six.moves.urllib.request import urlopen

    from sage.misc.sage_eval import sage_eval

    v = int(v)
    k = int(k)
    t = int(t)

    param = ("?v=%s&k=%s&t=%s"%(v,k,t))

    url = "http://www.ccrwest.org/cover/get_cover.php"+param
    if verbose:
        print "Looking up the bounds at %s" % url
    f = urlopen(url)
    s = f.read()
    f.close()

    if 'covering not in database' in s:   #not found
        str = "no (%d,%d,%d) covering design in database\n"%(v,k,t)
        raise ValueError(str)

    return sage_eval(s)
