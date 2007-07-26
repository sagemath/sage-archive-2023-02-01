
#*****************************************************************************
#           Copyright (C) 2007 Sean Howe and Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

class BruhatNode():
    """
    """
    def __init__(self):
        self.__permutation__ = None
        self.__degree__ = None
        self.__up_connected__ = [] # nodes that cover self
        self.__down_connected__ = [] # nodes that self covers

    def __get_copy__(self):
        copy=BruhatNode()
        copy.__permutation__=self.__permutation__
        return copy

class BruhatInterval():

    def __init__(self,height):
        self.__graph__ = None
        self.__nodes__ = []
        for i in range(0,height+1):
            self.__nodes__.append([])

    def create_graph(self):
        from sage.graphs import graph

        self.__graph__ = graph.DiGraph()
        for row in self.__nodes__:
                for node in row:
                    if node==None:
                        print len(row)
                    for node_above in node.__up_connected__:
                        self.__graph__.add_arc(node.__permutation__, node_above.__permutation__)

    def _sub_interval_add_below(self, sub_interval, to_add, height, added_by):
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
        sub_interval=BruhatInterval(height)

        x=None
        for x in self.__nodes__[height]:
            if end_perm==x.__permutation__:
                break
        self._sub_interval_add_below(sub_interval, x, height, None)

        return sub_interval

class BruhatSn(BruhatInterval):
    """
    """

    def __init__(self, n, height=None):
        if height is None:
            height = n*(n-1)/2

        BruhatInterval.__init__(self, height)

        first_perm=tuple(range(1,n+1))
        self._add_above(first_perm,0,height,None)

    def _swapitch(self, i, j, k):
        if k == i:
            return j
        elif k == j:
            return i
        else:
            return k

    def _add_above(self, permutation, height, max_height, added_by):
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




