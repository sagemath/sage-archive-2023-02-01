from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

cpdef list succ(x):
    cdef PermutationGroupElement u, u1, si
    cdef int n, first, i, j
    cdef list S, is_positive_root, successors

    S, n, is_positive_root, u, first = x
    successors = []
    for i in xrange(first):
        si = <PermutationGroupElement>(S[i])
        u1 = <PermutationGroupElement>(si._mul_(u))
        if test(is_positive_root,u1,i):
            successors.append((S, n, is_positive_root, u1, i))
    for i in xrange(first+1,n):
        if is_positive_root[u.perm[i]+1]:
            si = <PermutationGroupElement>(S[i])
            u1 = <PermutationGroupElement>(si._mul_(u))
            if test(is_positive_root,u1,i):
                successors.append((S, n, is_positive_root, u1, i))
    return successors

cdef bint test(list is_positive_root, PermutationGroupElement u1, int i):
    cdef int j

    for j in xrange(i):
        if not is_positive_root[u1.perm[j]+1]:
            return False
    return True

def search_forest_iterator(roots):
    cdef list stack

    stack = [iter(roots)]
    while stack:
        try:
            node = next(stack[-1])
        except StopIteration:
            stack.pop()
            continue

        yield node[3]
        stack.append( iter(succ(node)) )
