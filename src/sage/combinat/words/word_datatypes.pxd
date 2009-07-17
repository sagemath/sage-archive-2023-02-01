from sage.structure.sage_object cimport SageObject

# The generic class
cdef class WordDatatype(SageObject)

# finite WordDatatype
cdef class WordDatatype_list(WordDatatype)
cdef class WordDatatype_str(WordDatatype)
cdef class WordDatatype_tuple(WordDatatype)
