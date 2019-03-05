cdef extern from "hasse_diagram.h":
    #A pointer to the underlying C++ class of CombinatorialPolyhedron.
    ctypedef void* CombinatorialPolyhedron_ptr;
    
#    r"""
#    CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** facets_pointer, unsigned int nr_facets, unsigned int *len_facets, unsigned int nr_vertices, int is_unbounded)
    
#    Initalizes the C++ class of CombinatorialPolyhedron and returns a pointer to it.
    
#    INPUT:
#     - ``facets_pointer`` to a list of facets a lists of vertices/rays/lines.
#       The vertices/rays/lines need to be labeled 0,...,``nr_vertices - 1``.
#     - ``nr_facets`` the number of facets/the length of ``facets_pointer.
#     - ``len_facets`` the length of each face, i.e.
#       ``facetspointer[i]`` should be of length ``len_facets[i]``.
#     - ``is_unbounded`` needs to be 0, if the Polyhedron is unbounded,
#       otherwise 1 + ``nr_lines``.
#       ``is_unbounded = 0`` and ``is_unbounded = 1`` give the same result
#       for bounded Polyhedra, but ``is_unbounded = 0`` might be faster.
       
#    WARNING::
    
#        - ``nr_facets`` needs to be at least two.
#    """
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** facets_pointer, unsigned int nr_facets, unsigned int *len_facets, unsigned int nr_vertices, int is_unbounded)
    
#    r"""
#    CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** incidence_matrix, unsigned int nr_facets, unsigned int nr_vertices, int is_unbounded)
    
#    Initalizes the C++ class of CombinatorialPolyhedron and returns a pointer to it.
    
#    INPUT:
#     - An ``incidence_matrix`` given by a list of lenght ``nr_facets`` of
#       a list of ``nr_vertices`` of 0s and 1s.
#     - ``nr_facets`` the number of facets.
#     - ``nr_vertices`` the number of vertices/rays/lines.
#     - ``is_unbounded`` needs to be 0, if the Polyhedron is unbounded,
#       otherwise 1 + ``nr_lines``.
#       ``is_unbounded = 0`` and ``is_unbounded = 1`` give the same result
#       for bounded Polyhedra, but ``is_unbounded = 0`` might be faster.
    
#    WARNING::
    
#        - ``nr_facets`` needs to be at least two.
#    """
    cdef CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** incidence_matrix, unsigned int nr_facets, unsigned int nr_vertices, int is_unbounded)
    
#    r"""
#    delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)
    
#    Needs to be called to deallocate the memory used by the CombinatorialPolyhedron object in C++
#    """
    cdef void delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr)
    
#    r"""
#    unsigned int dimension(CombinatorialPolyhedron_ptr C)
    
#    Returns the dimension of the Polyhedron.
    
#    This is implicitly calculated and stored for most other functions.
#    """
    cdef unsigned int dimension(CombinatorialPolyhedron_ptr C)
    
#    r"""
#    void f_vector(CombinatorialPolyhedron_ptr C, unsigned long *vector)
    
#    Stores the f-vector of the Polyhedron in ``vector``. Vector needs to
#    be of length ``dimension + 2``
#    """
    cdef void f_vector(CombinatorialPolyhedron_ptr C, unsigned long *vector)
    
#    r"""
#    unsigned int ** edges(CombinatorialPolyhedron_ptr C)
    
#    Calculates the edges of the Polyhedron and returns a pointer to them.
    
#    The result ``edges`` will be a pointer of pointers to unsigned integers.
#    Each integer represents a vertex/ray/line in the order of the input at
#    initialization.
    
#    After getting the pointer, you can get the nr of edges by the corresponding entry in the f-vector.
    
#    ``edges[0]`` stores ``maxnumberedges`` edges.
#    (see ``get_maxnumberedges()``)
#    The first edge will have vertices ``edges[0][0]`` and ``edges[0][1]``,
#    the second edges ``edges[0][2]`` and ``edges[0][3]`` and so on.
#    The ``maxnumberedges + 1``st edge will have vertices
#    ``edges[1][0]`` and ``edges[1][1]``.
    
#    WARNING:
    
#        - ``edges`` will have at most length ``maxnumberedges``, so there will
#          be at most ``maxnumberedges^2`` edges stored.
        
#        - ``edges[1]`` is not initialized if the nr of edges is at most
#        ``maxnumberedges``.
#        - ``edges[i]`` is not initialized if the nr of edges is
#          too small.
    
#    NOTE: 
#        Might implicitly calculate the f-vector, but NOT vice versa.
        
#    """
    cdef unsigned int ** edges(CombinatorialPolyhedron_ptr C)
    
#    r"""
#    unsigned int ** ridges(CombinatorialPolyhedron_ptr C)
    
#    See ** edges(CombinatorialPolyhedron_ptr C).
#    """
    cdef unsigned int ** ridges(CombinatorialPolyhedron_ptr C)
    
#    r"""
#    unsigned long ** incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two, unsigned long * nr_incidences, unsigned int * twisted)
    
#    Calculates the incidences of ``dimension_one``-faces with ``dimension-two``-faces.
    
#    Return a pointer to ``incidences`` to the incedences similar to ``edges``.
    
#    INPUT:
    
#        - pointer to the C++-class of CombinatorialPolyhedron
#        - ``dimension_one`` and ``dimension_two`` the dimensions
#          of the faces of which the incidences will be given
#        - ``nr_incidences`` an ``unsigned long[1]`` in which the nr_of incidences will be stored
#        - ``twisted`` an ``unsigned int[1]``, which will indicate if the result is twisted/flipped
        
#    OUTPUT:
#        - ``incidences`` of type ``unsigned long[maxnumberincidences][maxnumberincidences]``.
#          ``incidences[0]`` stores ``maxnumberincidences`` incidenes.
#          (see ``get_maxnumberincidences()``)
#          incidences[1] stores another ``maxnumberincidences`` of incidences etc.
#        - ``twisted[0]`` indicates the output corresponding to an input 
#          where ``dimension_one`` and ``dimension_two`` are interchanged.
#        - ``nr_incidences[0]`` the number of incidences
        
#    The first incidences will have faces ``incidences[0][0]``
#    and ``incidences[0][1]``,
#    where the first is of dimension ``dimension_one`` and the second is
#    of length ``dimension_two``.
#    The faces are given as their index in the list of all faces as in 
#    ``get_faces``.
    
#    WARNING:
    
#        - ``incidences`` will have at most length ``maxnumberincidences``,
#          so there will be at most ``maxnumberincidences^2`` edges stored.
#        - ``incidences[1]`` is not initialized if the nr of incidences
#          at most ``maxnumberincidences``.
#        - ``incidences[i]`` is not initialized if the nr of incidences is
#          to small.
    
#    NOTE: 
#        Implicetely calculates a list of all faces of ``dimension_one``
#        and ``dimension_two``. AFTER calling ``incidences`` the order in
#        ``get_faces`` is the same as in ``face_iterator``.
    
#    """
    cdef unsigned long ** incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two, unsigned long * nr_incidences, unsigned int * twisted)
    
#    """
#    ``record_all_faces`` generates a list of all faces in C++ as to speed
#    up different calculations later.
#    """
    cdef void record_all_faces(CombinatorialPolyhedron_ptr C)
    
#    """
#    get_faces(CombinatorialPolyhedron_ptr C, int dimension, unsigned int facet_repr, unsigned int **faces_to_return, unsigned int *length_of_faces)
    
#    Gets a list of all faces.
    
#    INPUT:
    
#        - ``dimension`` of the faces to return.
#          Should be at least -1 and at most ``dimension.
#        - ``facet_repr`` -- 0 for the V-represenation of the faces.
#          1 for the H-representation.
#        - ``faces_to_return`` -- an array of type
#          ``unsigned int[nr_faces][nr_vertices]`` if ``facet_repr == 0``
#          and
#          ``unsigned int[nr_faces][nr_facets]`` if ``facet_repr == 1`]
#          NOTE:
#            - get the correct nr of faces by the entry in the f-vector
#            - nr_vertices is the nr of vertices/rays/lines, i.e. the 
#              length of the V-representation given on initialization
#            - If you really know what you are doing, you can initialize
#              ``faces_to_return[i]`` just long enough to make the longest
#              face fit.
#        - ``length_of_faces`` -- an array of typ ``unsigned int[nr_faces]
        
#    OUTPUT:
    
#        - faces_to_return will contain the faces as list of
#          vertices, if ``facet_repr == 0``, or
#          facets, if ``facet_repr == 1``.
#        - ``faces_to_return[i]`` will contain a face of length
#          ``length_of_faces[i]``.
#    """
    void get_faces(CombinatorialPolyhedron_ptr C, int dimension, unsigned int facet_repr, unsigned int **faces_to_return, unsigned int *length_of_faces)
    
#    """
#    face_iterator_init(CombinatorialPolyhedron_ptr C, int dimension, unsigned int vertex_repr, unsigned int facet_repr)
    
#    Will initalize a face_iterator.
    
#    INPUT:
    
#        - ``dimension`` of the faces to iterator over, if
#          ``dimension == -2`` then the iterator will use al faces
#        - ``vertex_repr`` needs to be set to 1 if vertex_repr is wanted,
#          0 otherwise
#        - ``facet_repr`` needs to beset to 1 if facet_repr is wanted,
#          both facet_repr and vertex_repr are possible

#    WARNING:
    
#        - There can be at most one working iterator around. The iterator
#          should be used up before calling different functions, especially,
#            - `faces`
#            - `vertices`
#            - `edges` (a second call is fine)
#            - `ridges` (a second call is fine)
#            - `f_vector` (a second call is fine)
#            - `record_all_faces`
#            - `incidences`
#            - `face_lattice`
#            - `flag`
#            - `k-simplicial` #not implemented yet
#            - `k-simple`    # not implemented yet
#    """
    cdef void face_iterator_init(CombinatorialPolyhedron_ptr C, int dimension, unsigned int vertex_repr, unsigned int facet_repr)

#    """
#    unsigned int face_iterator(CombinatorialPolyhedron_ptr C, unsigned int *Vface_to_return, unsigned int *Vlength, unsigned int *Hface_to_return, unsigned int *Hlength)
    
#    WARNING:
        
#        Do not call the face_iterator unless it has been initalized prior
#        by ``face_iterator_init``.
        
#        There can only be one face_iterator around.
    
#    INPUT:
    
#        - Pointer to CombinatorialPolyhedron in C.
#        - ``Vface_to_return`` of type ``unsigned int[nr_vertices]``.
#        - ``Vlength`` of type ``unsigned int[1]``
#        - ``Hface_to_return`` of type ``unsigned int[nr_facets]``.
#        - ``Hlength`` of type ``unsigned int[1]``.
        
#    OUTPUT:
    
#        - returns 0 if there was no next face, 1 otherwise
#        - if ``vertex_repr == 1`` (initalization)
#          ``Vface_to_return`` will contain ``Vlength[0]`` integers
#          corresponding to the vertices of the face
#        - if ``facet_repr == 1`` (initalization)
#          ``Hface_to_return`` will contain ``Hlength[0]`` integers
#          corresponding to the vertices of the face
#    """
    cdef unsigned int face_iterator(CombinatorialPolyhedron_ptr C, unsigned int *Vface_to_return, unsigned int *Vlength, unsigned int *Hface_to_return, unsigned int *Hlength)
    
    
    
#    """
#    get_flag(CombinatorialPolyhedron_ptr C, unsigned int *flagarray, unsigned int length)
    
#    Returns the number of flags of type ``flagarray``.
#    ``Flagarray`` needs to be a sorted array of length ``length``,
#    ``Flagarray contain integers from 0 up to ``dimension - 1``, where
#    ``dimensions`` is the dimension of the CombinatorialPolyhedron.
#    ``dimension = dimension(CombinatorialPolyhedron_ptr C)``
#    """
    cdef unsigned long get_flag(CombinatorialPolyhedron_ptr C, unsigned int *flagarray, unsigned int length)
    
    
#    """
#    ``get_maxnumberedges`` is a helper function for
#    ``edges(CombinatorialPolyhedron_ptr C)`` and
#    ``ridges(CombinatorialPolyhedron_ptr C)``
    
#    Edges and ridges are stored in an array of form
#    ``edges = unsigned int[maxnumberedges][maxnumberedges*2]``.
#    In order get the edges from the pointer, one needs this number.
#    ``edges`` does not store more than ``maxnumberedges^2`` edges.
#    """
    cdef unsigned long get_maxnumberedges()
    
#    """
#    ``get_maxnumberincidences`` is a helper function for
#    ``incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two, unsigned long * nr_incidences, unsigned int * twisted)``
    
#    Incidences are stored in an array of form
#    ``incidences = unsigned long[maxnumberincidences][maxnumberincidences*2]``.
#    In order get the incidences from the pointer, one needs this number.
#    ``incidences`` does not store more than ``maxnumberincidences^2`` incidences
#    """
    cdef unsigned long get_maxnumberincidences()
