
cpdef _flip_c(W, set positions, list extended_root_conf_indices, int i, side="both"):
    cdef int r, nr_ref, r_minus, j, k
    r = extended_root_conf_indices[i]
    nr_ref = W.nr_reflections()
    r_minus = (r + nr_ref) % (2*nr_ref) # get the negative root -r
    j = i
    for k in xrange(len(extended_root_conf_indices)):
        if extended_root_conf_indices[k] == r_minus and k not in positions and side != "positive":
            j = k
            break
        elif extended_root_conf_indices[k] == r and k not in positions and side != "negative":
            j = k
            break
    positions.remove(i)
    positions.add(j)
    if j != i:
        t = W.reflections()[ min(r,r_minus) ]
        for k in range(min(i,j)+1,max(i,j)+1):
            extended_root_conf_indices[k] = t( extended_root_conf_indices[k]+1 ) - 1
    return j

cpdef _construct_facets_c(list Q, w, int n=-1, int pos=0, int l=-1):
    r"""
        Returns the list of facets of the subword complex associated to the word Q and the element w in a Coxeter group W.
        """
    cdef int s
    cdef list X, Y
    if n == -1:
        n = len(Q)
    if l == -1:
        first = True
        l = w.length()
    else:
        first = False
    
    if l == 0:
        return [range(pos,n)]
    elif n < l+pos:
        return []
    
    s = Q[pos]
    X = _construct_facets_c(Q, w, n=n, pos=pos+1, l=l)
    for f in X:
        f.append(pos)
    
    if w.has_left_descent(s):
        Y = _construct_facets_c(Q, w.apply_simple_reflection_left(s), n=n, pos=pos+1, l=l-1)
        Y = X+Y
    else:
        Y = X
    if first:
        return sorted([ sorted(x) for x in Y ])
    else:
        return Y

