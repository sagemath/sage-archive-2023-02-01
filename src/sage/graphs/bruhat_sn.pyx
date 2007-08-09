"""
A module for dealing with the Bruhat ordering of the symmetric group.

AUTHORS:
    - Sean Howe, Emily Kirkman, Robert Miller (REU 2007, UW Seattle)
"""

#*****************************************************************************
#     Copyright (C) 2007 Sean Howe, Emily A. Kirkman and Robert L. Miller
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.graphs.graph import DiGraph
from sage.databases.database import SQLDatabase
from sage.dsage.dist_functions.dist_function import DistributedFunction
from sage.dsage.database.job import Job
import time
import os

###???
from sqlite3 import dbapi2 as sqlite
import re
from sage.databases.database import GenericSQLDatabase, SQLDatabase
###???

class BruhatSn(DiGraph):
    """
    Returns the Hasse diagram of the symmetric group Sn, with respect to
    Bruhat ordering.

    Lengths of permutations are stored, and an optional maximum length limits
    the construction to a neighborhood of the identity.

    """
    def __init__(self, int n, max_length=None, label=False):
        cdef int i

        if max_length is None:
            max_length = (n*(n-1))/2
        self.max_length = max_length
        self.min_length = 0
        self.length = max_length

        self.lengths = {}
        for i from 0 <= i <= max_length:
            self.lengths[i] = []

        self.identity = tuple(range(1,n+1))
        self.antidentity = tuple(range(n,0,-1))
        adj_dict = {}
        add_above(self, adj_dict, self.identity, None, 0, self.max_length, 0,
                  0, label)
        DiGraph.__init__(self, adj_dict)

    def plot(self, pos=None, layout=None, vertex_labels=True,
             edge_labels=False, vertex_size=200, graph_border=False,
             vertex_colors=None, partition=None, edge_colors=None,
             scaling_term=0.05, iterations=50,
             color_by_label=False, heights=None):
        """
        Overrides the normal plot function in GenericGraph, so most arguments
        are ignored.
        """
        G = self.to_undirected()
        return G.plot(heights=self.lengths, vertex_labels=False, vertex_size=5)

class BruhatIntervalSn(DiGraph):
    """
    The subinterval [start, end] = {v : start <= v <= end} of Sn under Bruhat
    ordering.

    """

    def __init__(self, start, end, check=True, label=False):
        cdef int a, b, n, i, j
        if check and not leq(start, end):
                raise TypeError("Must have start (%s) <= end (%s)"%(start, end))

        n = len(start)
        if check and n != len(end):
            raise TypeError("Start (%s) and end (%s) must have same length."\
                            %(start, end))

        a = length(start)
        b = length(end)
        self.max_length = max(a,b)
        self.min_length = min(a,b)
        self.length = self.max_length - self.min_length

        self.lengths = {}
        for i from 0 <= i <= self.max_length:
            self.lengths[i] = []

        self.identity = tuple(range(1,n+1))
        self.antidentity = tuple(range(n,0,-1))
        adj_dict = {}

        above_only = (end == self.antidentity)
        below_only = (start == self.identity)

        if above_only:
            add_above(self, adj_dict, start, None, a, self.max_length, 0, 0,
                      label)
        elif below_only:
            add_below(self, adj_dict, end, None, b, 0, 0, label)
        else:
            add_below(self, adj_dict, end, None, b, 0, 0, label)
            verts = vertices_above(start, adj_dict)
            for i from 0 <= i <= self.max_length:
                j = 0
                while j < len(self.lengths[i]):
                    u = self.lengths[i][j]
                    if u not in verts:
                        adj_dict.pop(u)
                        self.lengths[i].pop(j)
                    else:
                        j += 1
                        for v in adj_dict[u].keys():
                            if v not in verts:
                                adj_dict[u].pop(v)
        DiGraph.__init__(self, adj_dict)

        self.start = start
        self.end = end

    def is_self_dual(self):
        """
        Returns True iff interval is isomorphic to its dual, the same
        interval but with the order relation reversed.
        """
        return self.is_isomorphic(self.reverse())

    def plot(self, pos=None, layout=None, vertex_labels=False,
             vertex_size=200, graph_border=False,
             vertex_colors=None, partition=None, edge_colors=None,
             scaling_term=0.05, iterations=50,
             color_by_label=False, heights=None, **kwds):
        """
        Overrides the normal plot function in GenericGraph, so most arguments
        are ignored.
        """
        G = self.to_undirected()
        return G.plot(heights=self.lengths, vertex_labels=False,
                      vertex_size=5, **kwds)

def add_arc_bruhat(predecessor, successor, dict, i, j, label=False):
    if label: # know i < j and pred[i] < pred[j]
        # j - i == 1
        # <=> right by min(i,j)
        # <=> swap positions
        if j - i == 1:
            generator = [(i,'R')]
        else:
            generator = []
        # pred[j] - pred[i] == 1
        # <=> left by min(pred[i],pred[j])
        # <=> swap numbers
        if predecessor[j] - predecessor[i] == 1:
            generator.append((predecessor[i],'L'))
    else:
        generator = None
    dict[predecessor][successor] = generator

cdef void add_above(D, dict, current_perm, predecessor,
                    int current_length, int max_length, int a, int b,
                    label):
    """
    Recursive helper routine for contstructing BruhatSn.
    """
    cdef int i, j, k, n = len(current_perm)

    if current_perm in D.lengths[current_length]:
        add_arc_bruhat(predecessor, current_perm, dict, a, b, label)
        return
    # (else)
    D.lengths[current_length].append(current_perm)
    dict[current_perm] = {}
    if predecessor is not None:
        add_arc_bruhat(predecessor, current_perm, dict, a, b, label)
    else:
        dict[current_perm] = {}
    if current_length < max_length:
        for i from 0 <= i < n:
            for j from i < j < n:
                if current_perm[i] > current_perm[j]: continue
                safe = True
                for k from i < k < j:
                    if current_perm[i] < current_perm[k] and \
                       current_perm[k] < current_perm[j]:
                        safe = False
                        break
                if not safe: continue
                new_perm = transpose(current_perm, i, j, n)
                add_above(D, dict, new_perm, current_perm,
                          current_length+1, max_length, i, j, label)

cdef void add_below(D, dict, current_perm, successor,
                    int current_length, int a, int b, label):
    """
    Recursive helper routine for contstructing BruhatSn.
    """
    cdef int i, j, k, n = len(current_perm)

    if current_perm in D.lengths[current_length]:
        add_arc_bruhat(current_perm, successor, dict, a, b, label)
        return
    # (else)
    D.lengths[current_length].append(current_perm)
    dict[current_perm] = {}
    if successor is not None:
        add_arc_bruhat(current_perm, successor, dict, a, b, label)
    else:
        dict[current_perm] = {}
    for i from 0 <= i < n:
        for j from i < j < n:
            if current_perm[i] < current_perm[j]: continue
            safe = True
            for k from i < k < j:
                if current_perm[i] > current_perm[k] and \
                   current_perm[k] > current_perm[j]:
                    safe = False
                    break
            if not safe: continue
            new_perm = transpose(current_perm, i, j, n)
            add_below(D, dict, new_perm, current_perm,
                      current_length-1, i, j, label)

def vertices_above(v, dict, seen=None):
    """
    Returns the set of vertices u in dict (a dict of lists) such that u >= v.

    NOTE: not quite symmetric with vertices_below
    """
    if seen is None:
        seen = []
        vertices_above(v, dict, seen)
        return seen
    if v not in seen:
        seen.append(v)
        for u in dict[v]:
            vertices_above(u, dict, seen)

def vertices_below(v, D, seen=None):
    """
    Returns the set of vertices u in dict (a dict of lists) such that u >= v.

    NOTE: not quite symmetric with vertices_above
    """
    if seen is None:
        seen = []
        vertices_below(v, D, seen)
        return seen
    if v not in seen:
        seen.append(v)
        for u in D.predecessor_iterator(v):
            vertices_below(u, D, seen)

def transpose(perm, int i, int j, int n):
    """
    Swaps entries i and j in perm.
    """
    cdef int k
    new_perm = [0]*n
    for k from 0 <= k < n:
        if k == j:
            new_perm[k] = perm[i]
        elif k == i:
            new_perm[k] = perm[j]
        else:
            new_perm[k] = perm[k]
    return tuple(new_perm)

cdef int length(w):
    """
    Returns the length of permutation w.
    """
    cdef int i, n=len(w), L=0
    for i from 0 <= i < n:
        for j from i < j < n:
            if (w[i]>w[j]):
                L+=1
    return L

def leq(a, b):
    """
    Whether a <= b under Bruhat ordering.

    NOTE: False does not imply a > b (Bruhat ordering is not linear).
    """
    cdef int i, j, n, C, AA, BB
    n = len(a)
    A = [0 for _ in a]
    B = [0 for _ in a]
    for i from 0 <= i < n:
        AA = a[i]
        BB = b[i]
        for j from 0 <= j < i:
            if A[j] > AA:
                C = AA
                AA = A[j]
                A[j] = C
            if B[j] > BB:
                C = BB
                BB = B[j]
                B[j] = C
            if A[j] > B[j]: return False
        if AA > BB: return False
        A[i] = AA
        B[i] = BB
    return True

class BruhatDatabase(SQLDatabase):
    """
    Database for storing information about Bruhat intervals.

    """

    def __init__(self, filename):
        filename = os.path.realpath(filename)
        new = (not os.path.exists(filename))
        SQLDatabase.__init__(self, filename)
        if new:
            interval_classes_dict = {
            'class_label'   : {'sql':'TEXT',    'index':True,  'primary_key':True },
# label hash is INT not INTEGER, so that SQL does not automatically increment
            'length'        : {'sql':'INTEGER', 'index':True,  'primary_key':False},
            'num_verts'     : {'sql':'INTEGER', 'index':True,  'primary_key':False},
            'num_arcs'      : {'sql':'INTEGER', 'index':True,  'primary_key':False},
            'interval_type' : {'sql':'TEXT',    'index':True,  'primary_key':False},
            'self_dual'     : {'sql':'BOOLEAN', 'index':True,  'primary_key':False},
            'ranks'         : {'sql':'TEXT',    'index':True,  'primary_key':False},
            'rank_symmetric': {'sql':'BOOLEAN', 'index':True,  'primary_key':False}}

            intervals_dict = {
            'interval_label': {'sql':'INTEGER', 'index':True,  'primary_key':False},
            'start'         : {'sql':'TEXT',    'index':True,  'primary_key':False},
            'end'           : {'sql':'TEXT',    'index':True,  'primary_key':False}}

            permutations_dict = {
            'permutation'   : {'sql':'TEXT',    'index':True,  'primary_key':False},
            'disarray'      : {'sql':'INTEGER', 'index':False, 'primary_key':False},
            'skip'          : {'sql':'TEXT',    'index':False, 'primary_key':False},
            'interlock'     : {'sql':'INTEGER', 'index':False, 'primary_key':False},
            'block_not'     : {'sql':'TEXT',    'index':False, 'primary_key':False}}

            self.create_table('interval_classes', interval_classes_dict)
            self.create_table('intervals', intervals_dict)
            self.create_table('permutations', permutations_dict)

    def has_label_hash(self, label_hash):
        """
        Returns True if this hash has been seen before, False otherwise.

        """
        s = 'SELECT class_label FROM interval_classes WHERE class_label = '
        s += str(label_hash) # TODO
        b = self.__connection__.execute(s).fetchall()
        return len(b) > 0

    def has_permutation(self, perm):
        """
        Returns True if this permutation has been seen before.

        """
        s = 'SELECT permutation FROM permutations WHERE permutation = "'
        s += str(perm).replace(' ', '') + '"'
        b = self.__connection__.execute(s).fetchall()
        return len(b) > 0

    def permutation_hash(self, perm):
        """
        Returns the canonical label hash of perm.

        """
        s = 'SELECT interval_label FROM intervals WHERE start = '
        s += '"%s" AND '%(str(tuple(range(1,len(perm)+1))).replace(' ',''))
        s += 'end = "%s"'%(str(perm).replace(' ',''))
        b = self.__connection__.execute(s).fetchall()
        return b[0][0]

    def hashes_with_len_order_size(self, length, order, size):
        """
        Returns a list of canonical hashes that match the invariants.

        """
        s = 'SELECT class_label FROM interval_classes WHERE '
        s += 'length = %s AND num_verts = %s AND num_arcs = %s'\
             %(length, order, size)
        l = self.__connection__.execute(s).fetchall()
        return [a[0] for a in l]

    def id_hashes_with_len_order_size(self, length, order, size):
        """
        Returns a list of canonical hashes that match the invariants.

        """
        s = 'SELECT class_label FROM interval_classes WHERE '
        s += 'length = %s AND num_verts = %s AND num_arcs = %s'\
             %(length, order, size)
        s += ' AND interval_type = "id"'
        l = self.__connection__.execute(s).fetchall()
        return [a[0] for a in l]

    def commit_interval_class(self, label_hash, length, order, size,
                              interval_type, self_dual, ranks):
        """
        Records a row into the interval_classes table.

        INPUT:
            label_hash -- hash of the canonical label
            length -- int
            order -- int, number of vertices
            size -- int, number of arcs
            interval_type -- string
                Sn
                Identity
                Antidentity (we know b/c dual will be identity interval)
                Other
            self_dual -- bool
            ranks -- list of ranks: ranks[i] is the number of vertices of
                    length i

        """
        cdef int i = 0, j = len(ranks) - 1
        while i <= j and ranks[i] == 0:
            i += 1
        while i < j:
            if ranks[i] == ranks[j]:
                i += 1
                j -= 1
            else:
                i = j+2
        if i <= j+1:
            symmetric = 't'
        else:
            symmetric = 'f'
        if self_dual:
            self_dual = 't'
        else:
            self_dual = 'f'
        self.add_row( 'interval_classes',
                      (label_hash, length, order, size, interval_type,
                       self_dual, str(ranks).replace(' ', ''), symmetric),
                      ('class_label', 'length', 'num_verts', 'num_arcs',
                       'interval_type', 'self_dual', 'ranks', 'rank_symmetric') )

    def commit_interval(self, label_hash, start, end):
        """
        Records a row into the intervals table.

        INPUT:
            label_hash -- hash of the canonical label
            start -- tuple
            end -- tuple

        """
        self.add_row( 'intervals',
                      (label_hash, str(start).replace(' ', ''),
                       str(end).replace(' ', '')),
                      ('interval_label', 'start', 'end') )

    def commit_permutation(self, permutation, disarray, skip, interlock,
                           block_notation):
        """
        Records a row into the permutations table.

        INPUT:
            permutation -- tuple in map notation
            disarray -- int
            skip -- str
            interlock -- int
            block_notation -- str

        """
        self.add_row( 'permutations',
                      (str(permutation).replace(' ', ''), disarray, skip,
                       interlock, block_notation),
                      ('permutation', 'disarray', 'skip', 'interlock',
                       'block_not') )

    def publish_classes(self, directory):
        directory += '/'
        try:
            os.mkdir(directory)
        except:
            pass
        s = 'select class_label, length, num_verts, num_arcs, interval_type, self_dual, ranks, rank_symmetric from interval_classes'
        try:
            cur = self.__connection__.cursor()
            cur.execute(s)
            classes = cur.fetchall()
        except:
            raise RuntimeError('Failure to fetch query.')
        outlist = []
        i = 0
        html = []
        for interval_class in classes:
            output = ''
            i += 1
            label, length, order, size, type, dual, ranks, sym = interval_class
            s = 'select start, end from intervals where interval_label = "%s"'%label
            try:
                cur = self.__connection__.cursor()
                cur.execute(s)
                intervals = cur.fetchall()
            except:
                raise RuntimeError('Failure to fetch query.')
            start_rep = eval(intervals[0][0])
            end_rep = eval(intervals[0][1])
            G = BruhatIntervalSn(start_rep, end_rep)
            p = G.plot()
            p.save(directory + '%s.png'%i, figsize=[2,2])
            output += '\n    <tr>\n        <td bgcolor=white align=center rowspan=7>'
            output += '<img src="%s.png">\n        </td>'%i
            output += '\n        <td bgcolor=white align=right> Length of interval: \n        </td>'
            output += '\n        <td bgcolor=white align=left> %s \n        </td>\n    </tr>'%length
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Number of vertices: \n        </td>'
            output += '\n        <td bgcolor=white align=left> %s \n        </td>\n    </tr>'%order
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Number of edges: \n        </td>'
            output += '\n        <td bgcolor=white align=left> %s \n        </td>\n    </tr>'%size
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Type of interval: \n        </td>'
            output += '\n        <td bgcolor=white align=left> '
            if type == 'sn': output += 'Symmetric group'
            elif type == 'id': output += 'Identity interval'
            elif type == 'ai': output += 'Antidentity interval'
            else: output += 'Interior interval'
            output += '\n        </td>\n    </tr>'
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Self dual: \n        </td>'
            output += '\n        <td bgcolor=white align=left> '
            if dual == 't': output += 'True'
            else: output += 'False'
            output += ' \n        </td>\n    </tr>'
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Rank symmetric: \n        </td>'
            output += '\n        <td bgcolor=white align=left> '
            if sym == 't': output += 'True'
            else: output += 'False'
            output += ' \n        </td>\n    </tr>'
            output += '\n    <tr>'
            output += '\n        <td bgcolor=white align=right> Rank polynomial: \n        </td>'
            rank_poly = ''
            ranks = eval(ranks)
            minimum_rank = min([ii for ii in range(len(ranks)) if ranks[ii] != 0])
            for r in range(len(ranks)):
                if ranks[r] != 0:
                    if ranks[r] == 1:
                        coeff = ''
                    else:
                        coeff = str(ranks[r])
                    p = r - minimum_rank
                    if p == 0:
                        variable = ''
                    elif p == 1:
                        variable = 'q'
                    else:
                        variable = 'q<sup>%s</sup>'%p
                    expression = coeff + variable
                    if expression == '': expression = '1'
                    rank_poly += expression + ' + '
            rank_poly = rank_poly[:-3]
            output += '\n        <td bgcolor=white align=left> %s \n        </td>\n    </tr>'%rank_poly

            output += '\n    <tr>\n        <td colspan=3 bgcolor=white align=center>'
            output += 'Number of appearances in S<sub>n</sub>:\n        </td>\n    </tr>'
            appearances = {}
            minimum = -1
            for start, end in intervals:
                start = eval(start)
                end = eval(end)
                n = len(start)
                try:
                    appearances[n] += 1
                except:
                    appearances[n] = 1
                if n < minimum or minimum == -1:
                    minimum = n
                    start_rep = start
                    end_rep = end
            for n in appearances:
                output += '\n    <tr>\n        <td colspan=2 bgcolor=white align=right>'
                output += 'Appearances in S<sub>%s</sub>:\n        </td>'%n
                output += '\n        <td bgcolor=white align=left>%s'%appearances[n]
                output += '\n        </td>\n    </tr>'
            output += '\n    <tr>\n        <td colspan=3 bgcolor=white align=center>'

            output += '[ %s , %s ]\n        </td>\n    </tr>'%(start_rep, end_rep)

            output += '\n    <tr>\n        <td bgcolor=lightblue colspan=3 height=3>\n        </td>\n    </tr>'
            outlist.append([None,output])
        # sort outlist?
        output = '<html><table bgcolor=lightgrey cellpadding=3>\n'
        output += ''.join([out[1] for out in outlist])
        output += '\n</table></html>'
        f = file(directory + "index.html", 'w')
        f.write(output)
        f.close()

    def self_dual_permutations(self):
        pass

    def rank_symmetric_permutations(self):
        pass

##############################################################################

class DistributedBruhatIntervals(DistributedFunction):
    """
    DSAGE implementation of interval classification.

    INPUT:
        dsage -- running instance of dsage
        database -- a BruhatDatabase
        logfile -- file location & name

    """
    def __init__(self, dsage, database, name='BruhatIntervals',
                 logfile='log.txt', max_length=5, local_length=6):
        DistributedFunction.__init__(self, dsage)
        self.max_len = max_length
        self.loc_len = local_length
        self.cur_len = 1
        self.n = 2
        self.local_jobs = {}
        self.db = database
        self.name = name
        self.logfile = logfile
        self.log = open(logfile, 'a')
        self.log.write("\n\n%s - New DistributedBruhatIntervals Class\n\n"
                       %time.asctime())
        self.log.flush()

    def start(self):
        self.find_new_jobs()
        DistributedFunction.start(self)

    def restore(self, new_dsage):
        self.log = open(self.logfile, 'a')
        self.log.write('\n%s - Switching Over from Previous DSage Instance\n\n'
                       %time.asctime())
        self.log.flush()
        self.find_new_jobs()
        DistributedFunction.restore(self, new_dsage)

    def find_new_jobs(self):
        self.log.write('%s - Finding New Jobs\n'%time.asctime())
        self.log.flush()
        if self.cur_len < self.loc_len:
            # take care of locally, too small for dsage
            compute_locally_up_to_length(self.db, self.n, self.cur_len, self.loc_len)
        #else:
            # create jobs

def compute_locally_up_to_length(db, n, cur_len, max_len):
    while cur_len <= max_len:
        while n <= cur_len + 1:
            S = BruhatSn(n, max_len, label=True)
            for perm in S:
                if perm != S.identity:
                    compute_id_interval_locally(S, perm, db, cur_len)
            n += 1
        cur_len += 1
        # with this new cur_len, compute the minimum n needed to search.
        # if S_n is too small, there is no point. also n >= 2.
        n = 2
        while cur_len > (n*(n-1))/2:
            n += 1

def compute_locally_up_to_n(db, n, cur_len, max_n):
    max_len = (max_n*(max_n-1))/2
    while cur_len <= max_len:
        while n <= min(cur_len + 1,max_n):
            S = BruhatSn(n, max_len, label=True)
            for perm in S:
                if perm != S.identity:
                    compute_id_interval_locally(S, perm, db, cur_len)
            n += 1
        cur_len += 1
        # with this new cur_len, compute the minimum n needed to search.
        # if S_n is too small, there is no point. also n >= 2.
        n = 2
        while cur_len > (n*(n-1))/2:
            n += 1

def compute_id_interval_locally(S, perm, db, cur_len):
    """
    Makes Robert's eyes cross.

    """
    id_interval = None
    if not db.has_permutation(perm):
        id_interval = S.subgraph(vertices_below(perm, S))
        perm_result = analyze_permutation(id_interval, S.identity, perm)
        h_hash = perm_result[1]
        process_permutation(db, perm_result)
    else:
        h_hash = db.permutation_hash(perm)
    if perm in S.lengths[cur_len]:
        if id_interval is None:
            id_interval = S.subgraph(vertices_below(perm, S))
        num_verts = id_interval.order()
        num_arcs = id_interval.size()
        hashes = db.hashes_with_len_order_size(cur_len, num_verts, num_arcs)
        if h_hash not in hashes:
            id_hashes = db.id_hashes_with_len_order_size(cur_len, num_verts, num_arcs)
            eq_class_result = analyze_interval_class(S.identity, perm, h_hash, num_verts, num_arcs, cur_len, id_hashes)
            process_interval_class(db, eq_class_result)
        for v in id_interval:
            if v != S.identity and v != perm:
                interval_result = analyze_interval(v, perm)
                label_hash, length, order, size = process_interval(db, interval_result)
                hashes = db.hashes_with_len_order_size(length, order, size)
                if label_hash not in hashes:
                    id_hashes = db.id_hashes_with_len_order_size(length, order, size)
                    eq_class_result = analyze_interval_class(v, perm, label_hash, order, size, length, id_hashes)
                    process_interval_class(db, eq_class_result)

def analyze_interval_class(start, end, label_hash, order, size, length, id_hashes):
    cdef int i, n = len(start)
    H = BruhatIntervalSn(start, end)
    HR = H.reverse()
    HRC = HR.canonical_label()
    HRCH = HRC.dig6_string()
    self_dual = (HRCH == label_hash)
    ranks = []
    for i from 0 <= i <= H.max_length:
        try:
            ranks.append(len(H.lengths[i]))
        except:
            ranks.append(0)
    if tuple(range(1,n+1)) == start:
        if tuple(range(n,0,-1)) == end:
            type = 'sn'
        else:
            type = 'id'
    else:
        HR = H.reverse()
        if HRCH in id_hashes:
            type = 'ai'
        else:
            type = 'ot'
    return ('class', label_hash, length, order, size, type, self_dual, ranks)

def process_interval_class(db, result):
    _, label_hash, length, num_verts, num_arcs, type, self_dual, ranks = result
    # assert _ == 'class'
    db.commit_interval_class(label_hash, length, num_verts, num_arcs, type, self_dual, ranks)

def analyze_interval(start, end):
    G = BruhatIntervalSn(start, end) # TODO: check=False
    C = G.canonical_label()
    label_hash = C.dig6_string()
    return ('interval', label_hash, start, end, G.max_length - G.min_length, G.order(), G.size())

def process_interval(db, interval_result):
    _, label_hash, start, end, length, order, size = interval_result
    # assert _ == 'interval'
    db.commit_interval(label_hash, start, end)
    return label_hash, length, order, size

# alt: id = tuple(range(1,len(perm)+1))
#      id_interval = BruhatIntervalSn(id, perm, label=True)
def analyze_permutation(id_interval, start, end):
    cdef int i, disarray
    canonical_hash = id_interval.canonical_label().dig6_string()
    disarray = 0
    for i from 0 <= i < len(start):
        disarray += abs(i + 1 - end[i])
    skip = skip_from_permutation(end)
    interlock = interlock_from_permutation(end)
    block = blocks_from_permutation(end)
    block = str( block ).replace(' ', '')
    return ('permutation', canonical_hash, end, disarray, skip, interlock, block )

def process_permutation(db, perm_result):
    _, label_hash, perm, disarray, skip, interlock, block_not = perm_result
    # assert _ == 'permutation'
    db.commit_permutation(perm, disarray, skip, interlock, block_not)
    db.commit_interval(label_hash, tuple(range(1,len(perm)+1)), perm)

def blocks_from_permutation(perm):
    perm_blocks=[]
    cur_max = 1
    while cur_max <= len(perm):
        cur_block = []
        start_at = cur_max-1
        i = start_at
        while i < cur_max:
            cur_block.append(perm[i] - start_at)
            cur_max = max(cur_max, perm[i])
            i += 1
        if len(cur_block) != 1:
            perm_blocks.append(cur_block)
        cur_max += 1
    return perm_blocks

def num_cycles_in_block(block):
    n = len(block)
    orbits = [[0]]
    i = 0
    seen = [0]*n
    seen[0] = 1
    while i < n:
        dest = block[i]-1
        if seen[dest] == 1:
            i = 0
            while i < n and seen[i] == 1:
                i += 1
        else:
            for cell in orbits:
                if i in cell:
                    cell.append(dest)
                    i = dest
                    seen[i] = 1
                    break
        if i < n and seen[i] == 0:
            orbits.append([i])
    return len(orbits)

def skip_from_permutation(perm):
    perm_blocks = blocks_from_permutation(perm)
    blocks = []
    for x in perm_blocks:
        pattern = ''
        for i in range(len(x)):
            if x[i] == i+1:
                pattern += 'X'
            else:
                pattern += '_'
        split = pattern.split('X')
        gen_pattern = ''
        for i in range(1, len(split)):
            gen_pattern += 'X'
            if i != len(split)-1:
                gen_pattern += "%i"%(len(split[i]))
        if gen_pattern != '':
            blocks.append(gen_pattern)
    s = ''
    for x in blocks:
        s += '{' + x + '}'
    if s == '':
        s = '_'
    return s

def interlock_from_permutation(perm):
    perm_blocks = blocks_from_permutation(perm)
    interlock = 0
    for block in perm_blocks:
        interlock += num_cycles_in_block(block) - 1
    return interlock












#def vertices_below(v, D, seen=None):
#    """
#    Returns the set of vertices u in D such that u <= v.
#    """
#    if seen is None:
#        seen = []
#        vertices_below(v, D, seen)
#        return seen
#    if v not in seen:
#        seen.append(v)
#        for u in D.predecessor_iterator(v):
#            vertices_below(u, D, seen)
