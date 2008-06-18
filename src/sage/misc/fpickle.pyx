########################################################
# Function pickling
########################################################
import new, types, copy_reg, cPickle
#See python cookbook for more details
def code_ctor(*args):
    return new.code(*args)
def reduce_code(co):
    if co.co_freevars or co.co_cellvars:
        raise ValueError, "Cannot pickle code objects from closures"
    return code_ctor, (co.co_argcount, co.co_nlocals, co.co_stacksize,
                       co.co_flags, co.co_code, co.co_consts, co.co_names,
                       co.co_varnames, co.co_filename, co.co_name,
                       co.co_firstlineno, co.co_lnotab)

copy_reg.pickle(types.CodeType, reduce_code)

def pickle_function(func):
    """
    EXAMPLES:
        sage: def f(N): return N+1
        ...
        sage: g = pickle_function(f)
        sage: h = unpickle_function(g)
        sage: h(10)
        11
    """
    return cPickle.dumps(func.func_code)

def unpickle_function(pickled):
    """
    EXAMPLES:
        sage: def f(N,M): return N*M
        ...
        sage: unpickle_function(pickle_function(f))(3,5)
        15
    """
    recovered = cPickle.loads(pickled)
    ret = new.function(recovered, globals())
    return ret

