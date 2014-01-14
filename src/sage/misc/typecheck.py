################################################################################
#
# Method signature checking decorators
# Dmitry Dvoinikov <dmitry@targeted.org>
#
# Sample usage:
#
# class B(object):
#     @takes("B", basestring) # the third parameter (value) could have any type
#     def f(self, name, value):
#         return 0 # return type is not checked here
#
# class D(B):
#     @takes("D", basestring, B) # type D is unknown at this time and is provided by name
#     @returns(int)                  # whereas the third parameter could be explicitly typed as B
#     def f(self, name, value):
#         print name
#         return B.f(name, value) + 1 # has to return int
#
################################################################################

__all__ = ["takes", "InputParameterError", "returns", "ReturnValueError"]

################################################################################

def base_names(C):
    "Returns a tuple of names of base classes for a given class"

    return tuple([ x.__name__ for x in C.__mro__ ])

################################################################################

def takes(*types):
    "Method signature checking decorator"

    for t, i in zip(types, range(1, len(types) + 1)):
        if not isinstance(t, type) and not isinstance(t, str):
            raise TypeError("@takes decorator expects parameter %d to be a " \
                            "type or a string (type name), got %s" % \
                            (i, type(t).__name__))

    def takes_proxy(method):

        def takes_invocation_proxy(*args):

            if len(args) < len(types):
                raise InputParameterError("%s expects at least %d parameters, got %d" %
                                     (method.__name__, len(types), len(args)))

            for t, v, i in zip(types, args, range(1, len(args) + 1)):

                if isinstance(t, type) and not isinstance(v, t):
                    raise InputParameterError("%s expects parameter %d to have type %s, got %s" %
                                         (method.__name__, i, t.__name__, type(v).__name__))
                elif isinstance(t, str) and not t in base_names(v.__class__):
                    raise InputParameterError("%s expects parameter %d to have type %s, got %s" %
                                         (method.__name__, i, t, type(v).__name__))

            return method(*args)

        takes_invocation_proxy.__name__ = method.__name__ # just love Python for such tricks

        return takes_invocation_proxy

    return takes_proxy

class InputParameterError(TypeError): pass

################################################################################

def returns(sometype):
    "Return type checking decorator"

    if not isinstance(sometype, type) and not isinstance(sometype, str):
        raise TypeError("@returns decorator expected parameter to be a type " \
                        "or a string (type name), got %s" % type(sometype).__name__)

    def returns_proxy(method):

        def returns_invocation_proxy(*args):

            result = method(*args)

            if isinstance(sometype, type) and not isinstance(result, sometype):
                raise ReturnValueError("%s should have returned a value of type %s, got %s" %
                                       (method.__name__, sometype.__name__, type(result).__name__))
            elif isinstance(sometype, str) and not sometype in base_names(result.__class__):
                raise ReturnValueError("%s should have returned a value of type %s, got %s" %
                                       (method.__name__, sometype, type(result).__name__))
            else:
                return result

        returns_invocation_proxy.__name__ = method.__name__ # just love Python for such tricks

        return returns_invocation_proxy

    return returns_proxy

class ReturnValueError(TypeError): pass

################################################################################
