cimport cython

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef BraidOrbit(list word, list rels):
    cdef int l, rel_l, loop_ind, list_len
    cdef tuple left, right, test_word, new_word
    cdef list rel

    l = len(word)
    words      = set( [tuple(word)] )
    cdef list test_words = [ tuple(word) ]

    rels = rels + [ [b,a] for a,b in rels ]
    rels = [ [tuple(a), tuple(b), len(a)] for a,b in rels ]

    loop_ind = 0
    list_len = 1
    while loop_ind < list_len:
        test_word = <tuple> test_words[loop_ind]
        loop_ind += 1
        for rel in rels:
            left, right, rel_l = rel
            for i in range(l-rel_l+1):
                if pattern_match(test_word, i, left, rel_l):
                    new_word = test_word[:i]+right+test_word[i+rel_l:]
                    if new_word not in words:
                        words.add(new_word)
                        test_words.append(new_word)
                        list_len += 1
    return words

@cython.wraparound(False)
@cython.boundscheck(False)
cdef bint pattern_match(tuple L, int i, tuple X, int l):
    cdef int ind
    for ind in range(l):
        if not L[i+ind] == X[ind]:
            return False
    return True
