import graph
import bruhat
from sqlite3 import dbapi2 as sqlite
import os
import re
from sage.databases.database import GenericSQLDatabase, SQLDatabase
from sage.databases.db import DB_HOME
dblocation = DB_HOME + '/graphs/bruhat.db'

class BruhatDatabase(GenericSQLDatabase):
    """
    """

    def __init__(self):
        GenericSQLDatabase.__init__(self,dblocation)

    def harvest_data(self, start_height, end_height):
        """
        Returns a (mutable) copy of the database with rows added to all tables,
        reflecting the data obtained in the height range.

        Returns a SQLDatabase.

        TODO : This function is *unsafe* and should not remain in a released
        copy of sage.  (For development only).
        """
        db = self.copy()
        data_harvest(db,start_height,end_height)
        return db

    def commit_changes(self,db):
        """
        Assigns a copy of the database (with possibly additional data)
        to the location of BruhatDatabase.

        i.e.: Overwrites self with db in location.

        TODO : There is no reason for this function once harvest_data is
        removed.  Note:  This is here to break down the steps in case there is
        a failure in harvest_data
        """
        db.save(dblocation)

    def show(self, table_name, max_field_size=20, html_table=False, with_picture=None):
        if table_name is 'interval_classes':
            if with_picture is None:
                from sage.server.support import EMBEDDED_MODE
                with_picture = EMBEDDED_MODE
            if with_picture:
                from sage.plot.plot import plot
                from sage.structure.sage_object import load

                s = 'select picture, class_id, rep_interval, self_dual, kl_poly, atomic_number, skip_pattern, disarray, interlock, id_subintervals from interval_classes'

                try:
                    cur = self.__connection__.cursor()
                    cur.execute(s)
                    b = cur.fetchall()
                except:
                    raise RuntimeError('Failure to fetch query.')

                # get pictures
                for i in range (len(b)):
                    p = load(str(b[i][0]))
                    p.save('%s.png'%i, figsize=[2,2])

                print '<html><table bgcolor=lightgrey cellpadding=0><tr>'
                for i in range (len(b)):
                    print '<td bgcolor=white align=center rowspan="5"><img src="cell://%d.png"></td>'%i
                    print '<td bgcolor=white align=left> Class Id: %s </td>'%str(b[i][1])
                    print '<td bgcolor=white align=left> Atomic Number: %s </td>'%str(b[i][5])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left> Representative: %s </td>'%str(b[i][2])
                    print '<td bgcolor=white align=left> Skip Pattern: %s </td>'%str(b[i][6])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left> Self-Dual: %s </td>'%str(b[i][3])
                    print '<td bgcolor=white align=left> Disarray: %s </td>'%str(b[i][7])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left> K-L Poly: %s </td>'%str(b[i][4])
                    print '<td bgcolor=white align=left> Interlock: %s </td>'%str(b[i][8])
                    print '</tr><tr>'
                    print '<td bgcolor=white align=left colspan="2"> Subintervals from Identity: %s </td>'%str(b[i][9])
                    print '</tr><tr>'
                    if ( i != len(b)-1 ): print '<tr><td bgcolor=lightblue colspan="3" height="3"></td></tr>'
                print '</table></html>'

            else:
                GenericSQLDatabase.show(self, table_name, max_field_size, html_table)
        else:
            GenericSQLDatabase.show(self, table_name, max_field_size, html_table)

## MODULE

# ***********************************************************
# Compute Kaszdan-Lusztig Polynomials:
# ***********************************************************

def leq(a, b):
    i=0
    j=0
    n = len(a)
    A = [0 for _ in a]
    B = [0 for _ in a]
    for i in range(n):
        AA = a[i]
        BB = b[i]
        for j in range(i):
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

def l(w): # returns height of w
    word_length=len(w)
    height=0
    i=0
    for i in range(word_length):
        for j in range(i+1, word_length):
            if (w[i]>w[j]):
                height+=1
    return height

def get_first_descent(w): # returns first member of descent set of w
    word_length=len(w)
    for i in range(word_length-1):
        if (w[i]>w[i+1]):
            return [i, i+1]
    return []


def get_descent_set(w): # returns descent set of w
    word_length=len(w)
    descent_set=[]
    i=0
    for i in range(word_length-1):
        if (w[i]>w[i+1]):
            descent_set.append([i, i+1])
    return descent_set

def apply_transposition(w, s):
    new_w=w[:]
    new_w[s[0]]=w[s[1]]
    new_w[s[1]]=w[s[0]]

    return new_w

def get_immediately_above_set(w):
    set=[]
    length_w=len(w)
    for i in range(length_w):
        for j in range(i+1, length_w):
            if (w[i]<w[j]):
                works=true
                for k in range(i+1, j):
                    if (w[k]>w[i] and w[k]<w[j]):
                        works=false
                        break
                if (works==true):
                    above=w[:]
                    above[i]=w[j]
                    above[j]=w[i]
                    set.append(above)
    return set

def get_next_level(u, v):
    next_level=[]
    above_u=get_immediately_above_set(u)
    for x in above_u:
        if leq(x,v):
            next_level.append(x)
    return next_level

def get_set_inbetween(u,v):
    set=[]
    set.append([])
    set[0].append(u)

    start_height=l(u)
    end_height=l(v)
    cur_height=start_height
    while (cur_height<end_height):
        set.append([])
        index=cur_height-start_height
        for x in set[index]:
            above_x=get_next_level(x,v)
            for candidate in above_x:
                if candidate not in set[index+1]:
                    set[index+1].append(candidate)
        cur_height+=1
    return set

def mult_poly(a, b):
    c={}
    for x in a:
        for y in b:
            key=x+y
            if (c.has_key(key)):
                c[key]+=a[x]*b[y]
            else:
                c[key]=a[x]*b[y]
    return c

def add_poly(a,b):
    c={}
    for x in a:
        if (c.has_key(x)):
            c[x]+=a[x]
        else:
            c[x]=a[x]
    for x in b:
        if (c.has_key(x)):
            c[x]+=b[x]
        else:
            c[x]=b[x]
    return c

def display_poly(poly, var):
    str=""
    for x in poly:
        if (poly[x]!=0):
            str+="%.1f*%s^%.1f + " %(poly[x],var,x)
    str=str[:-3]
    print str

def get_poly_string(poly,var):
    str=""
    for x in poly:
        if (poly[x]!=0):
            str+="%.1f*%s^%.1f + " %(poly[x],var,x)
    str=str[:-3]
    return str

R_cache={}
def R_poly(u, v):
    hash=u.__str__()+";"+v.__str__()
    if (R_cache.has_key(hash)==true):
        return R_cache[hash]

    if (leq(u,v)!=True):
        R= {0:0}
    elif (u==v):
        R= {0:1}
    else:
        s=get_first_descent(v)
        new_u=apply_transposition(u, s)
        new_v=apply_transposition(v, s)

        du=get_descent_set(u)
        if (s in du):
            R= R_poly(new_u, new_v)
        else:
            R= add_poly(mult_poly({1:1},R_poly(new_u,new_v)), mult_poly({0:-1,1:1},R_poly(u,new_v)))
    R_cache[hash]=R
    return R

P_cache={}
def P_poly(x, w):
    hash=x.__str__()+";"+w.__str__()
    if (P_cache.has_key(hash)==true):
        return P_cache[hash]

    if (leq(x,w)==False):
        P= {0:0}
    elif (x==w):
        P= {0:1}
    else:
        l_w=l(w)
        l_x=l(x)

        sum={}
        inbetween=get_set_inbetween(x,w)[1:] # cull off the first node

        for level in inbetween:
            for y in level:
                 l_y=l(y)
                 term1= {0:(-1)^(l_y-l_x)}
                 #print term1
                 term2={((-l_x + 2*l_y -l_w)/2): 1}
                 term3=D_invol(R_poly(x,y))
                 term4=P_poly(y,w)
                 total=mult_poly(mult_poly(term1,term2),mult_poly(term3,term4))
                 sum=add_poly(sum,total)

        P=get_P_from_sum(sum,l_x,l_w)
       # print P
    P_cache[hash]=P
    return P

def D_invol(poly):
    D={}
    for x in poly:
        D[-x]=poly[x]
    return D

def get_P_from_sum(poly, l_x, l_w):
    P={}
    for x in poly:
        if (x<=0):
            new_key=x+((l_w-l_x)/2)
            P[new_key]=-poly[x]
    return P

def get_KL_poly(interval):
    bot = list(interval.__nodes__[0][0].__permutation__)
    top = list(interval.__nodes__[-1][0].__permutation__)
    return get_poly_string(P_poly(bot,top),'q')


# ***********************************************************
# Classify Interval and Permutation Classes
# ***********************************************************

def height(dg):
    d = {}
    for v in dg.vertices():
        h = height_of_vertex(v)
        try:
            d[h].append(v)
        except:
            d[h] = [v]
    return d

def height_of_vertex(v):
    h = 0
    n = len(v)
    for i in range(n):
        for j in range(i+1,n):
            if v[i] > v[j]:
                h += 1
    return h

def blocks_from_permutation(perm):
    perm_blocks=[]
    already_grouped=[]

    cur_max=1
    while cur_max<=len(perm):
        cur_block=[]
        start_at=cur_max-1
        i=start_at
        while i<cur_max:
            cur_block.append(perm[i]-start_at)
            cur_max=max(cur_max, perm[i])
            i+=1
        if (len(cur_block)!=1):
            perm_blocks.append(cur_block)
        cur_max+=1

    return perm_blocks

def cycles_from_permutation(perm):
    already_seen=[]
    cycles=[]
    for cur_index in range(0, len(perm)):
        if (perm[cur_index]-1 in already_seen):
            continue
        i=cur_index
        cur_cycle=[]
        while true:
            cur_cycle.append(i+1)
            already_seen.append(i)
            i=perm[i]-1
            if (i in already_seen):
                break
        if (len(cur_cycle)>1):
            cycles.append(cur_cycle)
    return cycles

def generators_from_permutation(perm):
    generators=""
    str=""
    cur_perm=[]
    for x in perm:
        cur_perm.append(x)
    while true:
        still_going=false
        for i in range(0, len(cur_perm)-1):
            if (cur_perm[i]>cur_perm[i+1]):
                swap_var=cur_perm[i]
                cur_perm[i]=cur_perm[i+1]
                cur_perm[i+1]=swap_var
                str="s_%i " %(i+1)
                generators=str+generators
                still_going=true
                break
        if (still_going==false):
            # We're at the base permutation
            break
    generators=generators[:len(generators)-1] #remove extra space at end
    return generators

def interval_classes_from_bruhat_graph(graph, max_height):
    class_list=[]

    for node in graph.__nodes__[max_height]:
        interval=graph.sub_interval_from_id(node.__permutation__,max_height)

        graph_already_classified=False
        interval_class=None
        for interval_class in class_list:
            if (interval_class.is_member(interval)):
                interval_class.add_occurrence(node.__permutation__)
                graph_already_classified=True
                break
        if (graph_already_classified==False):
            new_class=IntervalClass(interval, node.__permutation__)

            class_list.append(new_class)
    return class_list

def get_string_of_BruhatInterval(interval):
    perm = list(interval.__nodes__[-1][0].__permutation__)
    for i in range(len(perm)):
        if perm[len(perm)-1] == len(perm):
            perm.pop(-1)
        else:
            break
    return '[id, ' + str(tuple(perm)) + ']'

def get_interval_classes_at_height(height):
    n=2*height
    bruhat_sn=BruhatSn(n,height)
    ic=interval_classes_from_bruhat_graph(bruhat_sn, height)
    return ic

# ***********************************************************
# Database Modification:
# ***********************************************************

def add_interval_class_to_database(database, x, height, index, height_below): #height_below are the interval classes of the height below
    id=real("%i.%i" %(height,index+1))
    picture = dblocation + '/images/img%d_%d.sobj'%(height, index+1)
    # create picture and save it to that location:
    g = x.__representative_interval__.get_graph()
    G = g.to_undirected()
    p = G.plot(heights=height_unknown(g), vertex_labels=False, vertex_size=5)
    p.save(picture, figsize=[2,2])

    row = [picture, id, get_string_of_BruhatInterval(x.__representative_interval__), x.__self_dual__,get_KL_poly(x.__representative_interval__), x.__atomic_number__, x.__disarray__, x.__skip_pattern__.__str__(), x.__interlock__]

    subinterval_str=""
    for subinterval in x.__sub_intervals__:
        for i in range(len(height_below)):
            if (height_below[i].is_member(subinterval)):
                if not re.search('%d.%d,'%(height-1,i+1),subinterval_str):
                    subinterval_str+="%i.%i," %(height-1, i+1)
    subinterval_str=subinterval_str[:-1] # remove last comma
    row.append(subinterval_str)

    database.add_rows('interval_classes', [tuple(row)], \
                     ['picture','class_id','rep_interval','self_dual','kl_poly','atomic_number','disarray','skip_pattern','interlock','id_subintervals'])

    permRows=[]
    for pc in x.__perm_classes__:
        strBlocks=""
        for block in pc.__perm_blocks__:
            strBlocks+=block.get_perm_str()

        row=tuple([strBlocks, id])
        permRows.append(row)
    database.add_rows('permutation_classes', permRows, ['blocks','class_belongs_to'])

def data_harvest(database, start_height, end_height): #adds information to database and returns a list of interval classes at heights
    interval_classes_at_heights={}
    cur_height=start_height
    while cur_height <= end_height:
        interval_classes_at_heights[cur_height]=get_interval_classes_at_height(cur_height)

        for class_index in range(0, len(interval_classes_at_heights[cur_height])):
            if cur_height!=start_height:
                height_classes_below=interval_classes_at_heights[cur_height-1]
            else:
                height_classes_below=[]
            add_interval_class_to_database(database, interval_classes_at_heights[cur_height][class_index], cur_height, \
                                           class_index, height_classes_below)

        cur_height+=1

    return interval_classes_at_heights
