from sage.misc.all import tmp_filename, preparse

systems = {}
systems['PARI'] = ['sage.libs.pari', 'sage.interfaces.gp']
systems['Singular'] = ['sage.interfaces.singular', '_libsingular',
                       'sage.libs.singular']
systems['Maxima'] = ['sage.interfaces.maxima']
systems['GAP'] = ['sage.interfaces.gap']
systems['Magma'] = ['sage.interfaces.magma', 'sage.interfaces.magma_free']
systems['Axiom'] = ['sage.interfaces.axiom']
systems['ECM'] = ['sage.interfaces.ecm']
systems['scipy'] = ['scipy']
systems['numpy'] = ['numpy']
systems['Maple'] = ['sage.interfaces.maple']
systems['Mathematica'] = ['sage.interfaces.mathematica']
systems['MuPAD'] = ['sage.interfaces.mupad']
systems['Octave'] = ['sage.interfaces.octave']
systems['povray'] = ['sage.interfaces.povray']
systems['qsieve'] = ['sage.interfaces.qsieve']
systems['Macaulay2'] = ['sage.interfaces.macaulay2']
systems['mwrank'] = ['sage.interfaces.mwrank', 'sage.libs.mwrank']
systems['matlab'] = ['sage.interfaces.matlab']
systems['LiE'] = ['sage.interfaces.lie']
systems['Tachyon'] = ['sage.interfaces.tachyon']
systems['Frobby'] = ['sage.interfaces.frobby']
systems['gfan'] = ['sage.interfaces.gfan']
systems['R'] = ['sage.interfaces.r']
systems['KASH'] = ['sage.interfaces.kash']
systems['Linbox'] = ['sage.libs.linbox']
systems['Symmetrica'] = ['sage.libs.symmetrica']
systems['NTL'] = ['sage.libs.ntl', '_ntl',
                  'sage.rings.finite_field_ntl_gf2e']
systems['FLINT'] = ['_flint']
systems['GMP'] = ['sage.rings.integer.Integer']
systems['MPFR'] = ['sage.rings.real_mpfr',
                   'sage.rings.complex_number']
systems['MPFI'] = ['sage.rings.real_mpfi',
                   'sage.rings.complex_interval']
systems['M4RI'] = ['sage.matrix.matrix_mod2_dense']
systems['Givaro'] = ['sage.rings.finite_field_givaro']
systems['PolyBoRi'] = ['sage.rings.polynomial.pbori']

def get_systems(cmd):
    """
    Returns a list of the systems used in running the command
    cmd.  Note that the results can sometimes include systems
    that did not actually contribute to the computation. Due
    to caching and the inability to follow all C calls, it
    could miss some dependancies as well.

    INPUT:
       cmd -- a string to run

    EXAMPLES:
        sage: from sage.misc.citation import get_systems
        sage: get_systems('integrate(x^2, x)')
        ['Maxima']
        sage: R.<x,y,z> = QQ[]
        sage: I = R.ideal(x^2+y^2, z^2+y)
        sage: get_systems('I.primary_decomposition()')
        ['Singular']

    """
    import cProfile, pstats, re

    if not isinstance(cmd, basestring):
        raise TypeError, "command must be a string"

    cmd = preparse(cmd)

    #Run the command and get the stats
    filename = tmp_filename()
    cProfile.runctx(cmd, globals(), locals(), filename)
    stats = pstats.Stats(filename)

    #Strings is a list of method names and modules which get run
    strings = [a[0] + " " + a[2] for a in stats.stats.keys()]

    #Remove trivial functions
    bad_res = [re.compile(r'is_.*Element')]
    for bad_re in bad_res:
        i = 0
        while i < len(strings):
            if bad_re.findall(strings[i]):
                strings.pop(i)
            else:
                i += 1

    #Check to see which systems appear in the profiled run
    systems_used = []
    for system in systems:
        if any([(r in s) or (r.replace('.','/') in s) for r in systems[system] for s in strings]):
            systems_used.append(system)
    return systems_used
