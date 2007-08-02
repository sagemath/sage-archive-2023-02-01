"""
A module for dealing with Bruhat Intervals of the symmetric group.
Contains classes for BruhatNode, BruhatSn and BruhatInterval.

AUTHORS:
    - Sean Howe, Robert Miller, Emily Kirkman (2007-07-24): initial version

TODO:
 - all sub intervals (maybe own class)
"""
#*****************************************************************************
#           Copyright (C) 2007 Sean Howe, Emily A. Kirkman
#                                   and Robert L. Miller
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

class BruhatNode():
    """
    A node in a BruhatInterval.  Contains instance field variables to
    hold the permutation, degree and up and down neighbors of self.

    (Helper Class for BruhatInterval -- shouldn't get called on its own).
    """
    def __init__(self):
        self.__permutation__ = None
        self.__degree__ = None
        self.__up_connected__ = [] # nodes that cover self
        self.__down_connected__ = [] # nodes that self covers

    def __get_copy__(self):
        """
        Returns a deep copy of self without neighbor information.
        (For use in finding subintervals).
        """
        copy=BruhatNode()
        copy.__permutation__=self.__permutation__
        return copy

class BruhatInterval():
    """
    Class BruhatInterval.

    A class for constructing bruhat intervals on symmetric groups.
    Upon construction, a list structure is created to store the nodes
    of the bruhat interval at each height, and a graph variable is
    initialized to None, where the create_graph function will later
    store the DiGraph representing the interval.

    INPUT:
        end_height -- the greatest height of the interval
        start_height -- the least height of the interval (note that
            the identity is height=0)

    EXAMPLES:
        sage: b = BruhatSn(4)
        sage.: b.show_graph()
        sage: b_sub = b.sub_interval_from_id((4,2,3,1),5)
        sage.: b_sub.show_graph()
    """
    def __init__(self, end_height, start_height=None):
        if start_height is None:
            start_height = 0
        self.__graph__ = None
        self.__nodes__ = []
        for i in range(0,end_height+1):
            self.__nodes__.append([])

    def get_graph(self):
        """
        Returns the DiGraph representing the Bruhat Interval.

        EXAMPLES:
            sage: b = BruhatSn(4)
            sage.: (b.get_graph()).show()
            sage: b_sub = b.sub_interval_from_id((4,2,3,1),5)
            sage.: (b_sub.get_graph()).show()
        """
        return self.__graph__

    def _create_graph(self):
        """
        Creates a DiGraph representing the Bruhat Interval and stores
        it in self.__graph__
        """
        from sage.graphs import graph

        self.__graph__ = graph.DiGraph()
        for row in self.__nodes__:
                for node in row:
                    if node==None:
                        print len(row)
                    for node_above in node.__up_connected__:
                        self.__graph__.add_arc(node.__permutation__, node_above.__permutation__)

    def is_self_dual(self):
        """
        Returns true if the interval is self-dual, false if not
        """
        reverse=self.__graph__.reverse()
        return self.__graph__.is_isomorphic(reverse)

    def _sub_interval_add_below(self, sub_interval, to_add, height, added_by):
        """
        A recursive helper function for generating sub intervals.  (Should
        not be called directly).

        INPUT:
            sub_interval -- a BruhatInterval to add nodes to (recursively)
            to_add -- the BruhatNode to add next, which will be copied and
                added to the current sub_interval, unless it is already in
                the sub_interval; in which case a connection will be added
                between the copy and added_by
            height -- the current height (i.e.: the height of to_add)
            added_by -- a BruhatNode in sub_interval that is immediately
                above to_add
        """
        alreadyAdded=False
        x=None
        for x in sub_interval.__nodes__[height]:
            if x.__permutation__==to_add.__permutation__:
                alreadyAdded=True
                break
        if (alreadyAdded==True):
            if (added_by!=None):
                x.__up_connected__.append(added_by)
                added_by.__down_connected__.append(x)
            return

        adding=to_add.__get_copy__()
        sub_interval.__nodes__[height].append(adding)
        if (added_by!=None):
            adding.__up_connected__.append(added_by)
            added_by.__down_connected__.append(adding)

        for child in to_add.__down_connected__:
            self._sub_interval_add_below(sub_interval,child,height-1,adding)

    def sub_interval_from_id(self, end_perm, height):
        """
        Returns a BruhatInterval instance for the subinterval of self
        that starts at the identity and ends at end_perm.

        INPUT:
            end_perm -- the permutation at maximal height in the sub interval
            height -- the height of end_perm

        NOTE:
            If it is known that you will start from the identity, using this
            method is much faster than the general sub_interval function.

        EXAMPLES:
            sage: b = BruhatSn(4)
            sage.: b.show_graph()
            sage: b_sub = b.sub_interval_from_id((4,2,3,1),5)
            sage.: b_sub.show_graph()
        """
        sub_interval=BruhatInterval(height)

        x=None
        for x in self.__nodes__[height]:
            if end_perm==x.__permutation__:
                break
        self._sub_interval_add_below(sub_interval, x, height, None)

        sub_interval._create_graph()
        return sub_interval

    def is_isomorphic(self, other):
        """
        The comparison operator for two instances of BruhatInterval.
        Returns True if the two intervals are isomorphic, and False
        otherwise.

        Relies on NICE algorithm isomorphism tester.

        INPUT:
            other -- an instance of BruhatInterval to test for
                isomorphism with self

        EXAMPLES:
            sage: b = BruhatSn(4)
            sage.: b.show_graph()
            sage: b_sub = b.sub_interval_from_id((4,2,3,1),5)
            sage.: b.show_graph()
            sage: b.is_isomorphic(b_sub)
            False
            sage: a1 = b.sub_interval_from_id((2,3,1,4),2)
            sage: a2 = b.sub_interval_from_id((3,1,2,4),2)
            sage: a1.is_isomorphic(a2)
            True
        """
        g = self.__graph__
        h = other.__graph__

        # TODO : When NICE is updated, remove the following 4 lines:
        if g.size() != h.size():
            return False
        if g.order() != h.order():
            return False

        return g.is_isomorphic(h)

    def show_graph(self, labels=False):
        """
        Shows the graph nicely (i.e.: with nodes of same height drawn at
        same height).

        INPUT:
            labels -- whether or not to display node labels (permutations)

        EXAMPLES:
            sage: b = BruhatSn(4)
            sage.: b.show_graph()
            sage: b_sub = b.sub_interval_from_id((4,2,3,1),5)
            sage.: b_sub.show_graph()
        """
        heights=height_unknown(self.__graph__)
        self.__graph__.show(heights=heights, vertex_size=10, vertex_labels=labels)

class BruhatSn(BruhatInterval):
    """
    Class BruhatSn.

    Constructs the Bruhat Interval (or a portion of it up to user-defined
    height) from the symmetric group.

    INPUT:
        n -- As in S_n !
        height -- the height to build the interval up to.  (If left as
            default, then BruhatSn will construct the entire group).

    EXAMPLES:
        sage: glist = []
        sage: for i in range(2,6):
        ...       b = BruhatSn(i)
        ...       glist.append(b.get_graph())
        ...
        sage.: graphs_list.show_graphs(glist, layout='spring')
    """
    def __init__(self, n, height=None):
        if height is None:
            height = n*(n-1)/2

        BruhatInterval.__init__(self, height)

        first_perm=tuple(range(1,n+1))
        self._add_above(first_perm,0,height,None)
        self._create_graph()

    def _swapitch(self, i, j, k):
        """
        TODO : Robert - What does this do????
        """
        if k == i:
            return j
        elif k == j:
            return i
        else:
            return k

    def _add_above(self, permutation, height, max_height, added_by):
        """
        Recursive helper function for adding nodes above current place
        in interval.

        INPUT:
            permutation -- the current permutation
            height -- the current height
            max_height -- the stopping height
            added_by -- the node below
        """
        permutation_length = len(permutation)
        already_added = False

        for node in self.__nodes__[height]:
            if node.__permutation__ == permutation:
                already_added = True
                if (added_by != None):
                    node.__down_connected__.append(added_by)
                    added_by.__up_connected__.append(node)
                break

        if already_added == False:
            new_node = BruhatNode()
            new_node.__permutation__ = permutation
            new_node.__height__ = height
            if (added_by != None):
                new_node.__down_connected__.append(added_by)
                added_by.__up_connected__.append(new_node)
            self.__nodes__[height].append(new_node)

            if height < max_height:
                for i in range(0, permutation_length):
                    for j in range(i+1, permutation_length):
                        if permutation[i] >= permutation[j]: continue
                        safe = True
                        for k in range(i+1, j):
                            if (permutation[i] < permutation[k]) and (permutation[k] < permutation[j]):
                                safe = False
                                break
                        if not safe: continue

                        w = tuple([permutation[self._swapitch(i,j,k)] for k in range(0,permutation_length)])
                        self._add_above(w, height+1, max_height, new_node)

class IntervalClass:
    def __init__(self, interval, first_occurrence):
        self.__representative_interval__=interval

        #check if self-dual
        self.__self_dual__=self.__representative_interval__.is_self_dual()

        # find sub intervals
        self.__sub_intervals__=[]
        height_of_children=len (interval.__nodes__) -2
        for child in interval.__nodes__[height_of_children]:
            self.__sub_intervals__.append(interval.sub_interval_from_id(child.__permutation__, height_of_children))

        # add first occurrence, get invariants
        self.__perm_classes__=[]
        self.add_occurrence(first_occurrence)
        self.__disarray__=self.__perm_classes__[0].get_disarray()
        self.__atomic_number__=self.__perm_classes__[0].get_atomic_number()
        self.__interlock__=self.__perm_classes__[0].get_interlock()
        self.__skip_pattern__=self.__perm_classes__[0].get_skip_pattern()

    def is_member(self,other_interval):
        return other_interval.is_isomorphic(self.__representative_interval__)

    def add_occurrence(self,end_perm):
        class_exists=False
        for pc in self.__perm_classes__:
            if (pc.is_member(end_perm)):
                pc.add_occurrence(end_perm)
                class_exists=True
                break

        if (class_exists==False):
            newClass=PermutationClass(end_perm)
            self.__perm_classes__.append(newClass)


    def __str__(self):
        str=""
        str+="Self-dual: "
        if (self.__self_dual__==True):
            str+="True\n"
        else:
            str+="False\n"
        str+="Disarray:%i" %(self.__disarray__) +"\n"
        str+="Atomic Number:%i" %(self.__atomic_number__) +"\n"
        str+="Interlock:%i" %(self.__interlock__) +"\n"
        str+="Skip Pattern:" + (self.__skip_pattern__.__str__()) + "\n"
        str+="Permutation Classes:\n"
        for pc in self.__perm_classes__:
            str+= "\n"+pc.__str__()
        return str

class IntervalOccurrence:
    def __init__(self, perm):
        self.__end_perm__=perm
        self.__perm_length__=len(perm)

class SkipPattern:
    def __init__(self,perm_blocks):
        self.__blocks__=[]

        for x in perm_blocks:
            pattern=""
            for i in range(len(x.__perm__)):
                if (x.__perm__[i])==i+1:
                    pattern+="X"
                else:
                    pattern+="_"
            split=pattern.split('X')
            gen_pattern=""
            for i in range(1, len(split)):
                gen_pattern+="X"
                if (i!=len(split)-1):
                    gen_pattern+="%i" %(len(split[i]))

            if gen_pattern!="":
                self.__blocks__.append(gen_pattern)

    def is_equal(self, other):
        if len(other.__blocks__)!=len(self.__blocks__):
            return False

        used=[]
        for i in range(len(other.__blocks__)):
            used[i]=False

        for i in range(len(self.__blocks__)):
            matched=False
            for j in range(len(other.__blocks__)):
                if used[j]==True:
                    continue
                elif self.__blocks__[i]==other.__blocks__[j]:
                    used[j]=True
                    matched=True
                    break
                if matched==False:
                    return False
        return True

    def __str__(self):
        str=""
        for x in self.__blocks__:
            str+="{"+x+"}"
        if str=="":
            return "_"
        return str

class PermutationBlock:
    def __init__(self, perm_block):
        self.__perm__=perm_block[:]
        self.__cycles__=cycles_from_permutation(self.__perm__)
        self.__generators__=generators_from_permutation(self.__perm__)

    def is_equal(self,other_perm):
        return self.__perm__==other_perm

    def get_skip_pattern(self):
        return SkipPattern(self.__perm__)

    def get_cycle_str(self):
        str=""
        for cycle in self.__cycles__:
            str+="("
            for x in cycle:
                str+="%i " %(x)
            str=str[:len(str)-1]+")"
        return str

    def get_perm_str(self):
        str="["
        for x in self.__perm__:
            str+="%i " %(x)
        str=str[:len(str)-1]+"]"
        return str

    def get_generator_str(self):
        return self.__generators__

class PermutationClass:
    def __init__(self, perm):
        self.__perm_blocks__=[]
        blocks=blocks_from_permutation(perm)
        for block in blocks:
            new_pb=PermutationBlock(block)
            self.__perm_blocks__.append(new_pb)

        self.__occurrences__=[]
        self.add_occurrence(perm)

    def is_member(self,other_perm):
        blocks=blocks_from_permutation(other_perm)
        if (len(blocks)!=len(self.__perm_blocks__)):
            return False
        used=[]
        for block in blocks:
            used.append(False)

        for i in range(0,len(self.__perm_blocks__)):
            matched=False
            for j in range(0,len(blocks)):
                if (used[j]==False):
                    if(self.__perm_blocks__[i].is_equal(blocks[j])):
                        matched=True
                        used[j]=True
                        break
            if (matched==False):
                return False

        return True

    def get_disarray(self):
        disarray=0
        for block in self.__perm_blocks__:
            perm=block.__perm__
            for x in range(1,len(perm)+1):
                for i in range(0, len(perm)):
                    if (x==perm[i]):
                        disarray+=abs((i+1)-x)
                        break
        return disarray

    def get_atomic_number(self):
        atomic=0
        for block in self.__perm_blocks__:
            block_atomic=0
            for cycle in block.__cycles__:
                block_atomic+=len(cycle)
            atomic+=block_atomic-1
        return atomic

    def get_interlock(self):
        interlock=0
        for block in self.__perm_blocks__:
            interlock+=len(block.__cycles__)-1
        return interlock

    def add_occurrence(self,perm):
        occurrence=IntervalOccurrence(perm)
        self.__occurrences__.append(occurrence)

    def get_skip_pattern(self):
        pattern=SkipPattern(self.__perm_blocks__)

        return pattern

    def __str__(self):
        str=""
        for block in self.__perm_blocks__:
            str+=block.get_perm_str()
        str+=" = "
        for block in self.__perm_blocks__:
            str+="{"
            str+=block.get_cycle_str()
            str+="}"
        str+=" = "
        for block in self.__perm_blocks__:
            str+="{"+block.get_generator_str()+"}"
        return str

def height_unknown(G):
    """
    Helper function to determine height for plotting.
    """
    h = {}
    for v in G:
        if G.in_degree(v) == 0:
            break
    i = 0
    h = recurse_unknown(G, i, v, h)
    return h

def recurse_unknown(G, i, v, h, seen=None):
    """
    Helper function to determine height for plotting.
    """
    if seen is None: seen = []
    seen.append(v)
    try:
        h[i].append(v)
    except:
        h[i] = [v]
    for a, b, _ in G.outgoing_arc_iterator(v):
        if b not in seen:
            h = recurse_unknown(G, i+1, b, h, seen)
    return h


