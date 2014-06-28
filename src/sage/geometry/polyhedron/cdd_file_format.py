"""
Generate cdd ``.ext`` / ``.ine`` file format
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from misc import _set_to_None_if_empty, _common_length_of, _to_space_separated_string



#########################################################################
def cdd_Vrepresentation(cdd_type, vertices, rays, lines):
    r"""
    Return a string containing the V-representation in cddlib's ext format.

    NOTE:

    If there is no vertex given, then the origin will be implicitly
    added. You cannot write the empty V-representation (which cdd
    would refuse to process).

    EXAMPLES::

        sage: from sage.geometry.polyhedron.cdd_file_format import cdd_Vrepresentation
        sage: print cdd_Vrepresentation('rational', [[0,0]], [[1,0]], [[0,1]])
        V-representation
        linearity 1 1
        begin
          3 3 rational
          0 0 1
          0 1 0
          1 0 0
        end
    """
    vertices = _set_to_None_if_empty(vertices)
    rays     = _set_to_None_if_empty(rays)
    lines    = _set_to_None_if_empty(lines)

    num, ambient_dim = _common_length_of(vertices, rays, lines)

    # cdd implicitly assumes that the origin is a vertex if none is given
    if vertices is None:
        vertices = [[0]*ambient_dim]
        num += 1

    s = 'V-representation\n'
    if lines is not None:
        n = len(lines)
        s += "linearity " + repr(n) + ' '
        s += _to_space_separated_string(range(1,n+1)) + '\n'
    s += 'begin\n'
    s += ' ' + repr(num) + ' ' + repr(ambient_dim+1) + ' ' + cdd_type + '\n'
    if lines is not None:
        for l in lines:
            s += ' 0 ' + _to_space_separated_string(l) + '\n'
    if rays is not None:
        for r in rays:
            s += ' 0 ' + _to_space_separated_string(r) + '\n'
    if vertices is not None:
        for v in vertices:
            s += ' 1 ' + _to_space_separated_string(v) + '\n'
    s += 'end\n'
    return s

#########################################################################
def cdd_Hrepresentation(cdd_type, ieqs, eqns):
    r"""
    Return a string containing the H-representation in cddlib's ine format.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.cdd_file_format import cdd_Hrepresentation
        sage: cdd_Hrepresentation('rational', None, [[0,1]])
        'H-representation\nlinearity 1 1\nbegin\n 1 2 rational\n 0 1\nend\n'
    """
    ieqs = _set_to_None_if_empty(ieqs)
    eqns  = _set_to_None_if_empty(eqns)

    num, ambient_dim = _common_length_of(ieqs, eqns)
    ambient_dim -= 1

    s = 'H-representation\n'
    if eqns is not None:
        assert len(eqns)>0
        n = len(eqns)
        s += "linearity " + repr(n) + ' '
        s += _to_space_separated_string(range(1,n+1)) + '\n'
    s += 'begin\n'
    s += ' ' + repr(num) + ' ' + repr(ambient_dim+1) + ' ' + cdd_type + '\n'
    if eqns is not None:
        for e in eqns:
            s += ' ' + _to_space_separated_string(e) + '\n'
    if ieqs is not None:
        for i in ieqs:
            s += ' ' + _to_space_separated_string(i) + '\n'
    s += 'end\n'
    return s
