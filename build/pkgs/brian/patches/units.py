# ----------------------------------------------------------------------------------
# Copyright ENS, INRIA, CNRS
# Contributors: Romain Brette (brette@di.ens.fr) and Dan Goodman (goodman@di.ens.fr)
#
# Brian is a computer program whose purpose is to simulate models
# of biological neural networks.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# ----------------------------------------------------------------------------------
#
######### PHYSICAL UNIT CLASSES ###################
#------------------------------------ Dan Goodman -
#todo: More additional derived units (for the lookup table)

"""Defines physical units and quantities

The standard way to use this class is as follows:

V = 3 * volt
I = 2 * amp
R=V/I
print R

will return

1.5 ohm

The following fundamental units are defined:

metre, kilogram, second, amp, kelvin, mole, candle

And these additional basic units:

radian, steradian, hertz, newton, pascal, joule, watt,
coulomb, volt, farad, ohm, siemens, weber, tesla, henry,
celsius, lumen, lux, becquerel, gray, sievert, katal,
gram, gramme

Additionally, it includes all scaled versions of these
units using the standard SI prefixes (see the documentation
for the Unit class for more details), e.g. uamp,
mmetre, etc. It also includes the second and third powers
of each of these units, e.g. mvolt2 = mvolt*mvolt,
metre3 = metre**3, etc.

The module also defines these classes:

-- Dimension
        Stores the physical dimensions (length, mass, etc.)
-- DimensionMismatchError
        Exception raised if you try to add inconsistent units,
        etc.
-- Quantity
        The class of a value with a unit
-- Unit
        The class of the defined units like mvolt, etc.
-- UnitRegistry
        Stores 'known' units for printing

These functions:

-- get_dimensions(x)
        Returns the dimensions of a quantity or number x
-- have_same_dimensions(x,y)
        Tests if x and y have the same dimensions
-- is_dimensionless(x)
        Tests if x is dimensionless
-- display_in_unit(x,u)
        Displays quantity x in unit u
-- register_new_unit(u)
        Add a new unit u to the list of 'known' units for
        printing purposes
-- get_unit(x)
        Returns the fundamental unit of value x if one is known, or
        simply the value 1 with dimensions of x if none is known

And this decorator for function argument checking:

-- check_units(...)

If you want to use shorter named units, import the stdunits
module, which defines things like mV for mvolt, etc. They
are not included by default in the units module because of
the potential for variable name clashes.
"""

__origall__ = ['Dimension', 'Scale', 'DimensionMismatchError',
           'get_dimensions', 'is_dimensionless', 'have_same_dimensions',
           'display_in_unit', 'Quantity', 'Unit', 'register_new_unit',
           'check_units', 'is_scalar_type', 'get_unit', 'get_unit_fast',
           'scalar_representation']

__all__ = __origall__ + []

from brian_unit_prefs import bup
from operator import isNumberType, isSequenceType
from itertools import izip
import math, numpy
from utils.approximatecomparisons import *
import types
from functools import *
import sys
# Note that the decorator module below is used to provide signature preserving
# decorators, but it has the unfortunate side effect of messing up the tracebacks
# because it uses eval, so we only use it when we want to generate documentation,
# i.e. if 'sphinx' or 'docutils' or 'epydoc' are loaded.
try:
    import decorator
    use_decorator = 'sphinx' in sys.modules or 'docutils' in sys.modules or 'epydoc' in sys.modules
except:
    use_decorator = False

# SI dimensions (see table at end of file) and various descriptions,
# each description maps to an index i, and the power of each dimension
# is stored in the variable dims[i]
_di = { "Length":0, "length": 0, "metre":0, "metres":0, "metre": 0, "metres":0, "metre":0, "metres":0, "metre": 0, "metres":0, "m": 0, \
       "Mass":1, "mass": 1, "kilogram":1, "kilograms":1, "kilogram": 1, "kilograms":1, "kg": 1, \
       "Time":2, "time": 2, "second":2, "seconds":2, "second": 2, "seconds":2, "s": 2, \
       "Electric Current":3, "Electric Current":3, "electric current": 3, "Current":3, "current":3, "ampere":3, "amperes":3, "ampere": 3, "amperes":3, "A": 3, \
       "Temperature":4, "temperature": 4, "kelvin":4, "kelvins":4, "kelvin": 4, "kelvins":4, "K": 4, \
       "Quantity of Substance":5, "Quantity of substance": 5, "quantity of substance": 5, "Substance":5, "substance":5, "mole":5, "moles":5, "mole": 5, "moles":5, "mol": 5, \
       "Luminosity":6, "luminosity": 6, "candle":6, "candles":6, "candle": 6, "candles":6, "cd": 6 }
_ilabel = ["m", "kg", "s", "A", "K", "mol", "cd"]
# SI unit _prefixes, see table at end of file
_siprefixes = {"y":1e-24, "z":1e-21, "a":1e-18, "f":1e-15, "p":1e-12, "n":1e-9, "u":1e-6, "m":1e-3, "c":1e-2, "d":1e-1, \
            "":1, \
            "da":1e1, "h":1e2, "k":1e3, "M":1e6, "G":1e9, "T":1e12, "P":1e15, "E":1e18, "Z":1e21, "Y":1e24}


class Dimension(object):
    '''Stores the indices of the 7 basic SI unit dimension (length, mass, etc.)

    Provides a subset of arithmetic operations appropriate to dimensions:
    multiplication, division and powers, and equality testing.

    Methods:

    is_dimensionless() returns Boolean value

    Notes:

    Most users shouldn't use this class directly, but instead write things
    like:

    x = 3 * mvolt, etc.
    '''
    __slots__ = ["_dims"]
    #### INITIALISATION ####
    def __init__(self, *args, **keywords):
        """Initialise Dimension object with a vector or keywords

        Call as Dimension(list), Dimension(keywords) or Dimension(dim)

        list -- a list with the indices of the 7 elements of an SI dimension
        keywords -- a sequence of keyword=value pairs where the keywords are
          the names of the SI dimensions, or the standard unit
        dim -- a dimension object to copy

        Examples:

        The following are all definitions of the dimensions of force

        Dimension(length=1, mass=1, time=-2)
        Dimension(m=1, kg=1, s=-2)
        Dimension([1,1,-2,0,0,0,0])

        The 7 units are (in order):

        Length, Mass, Time, Electric Current, Temperature,
        Quantity of Substance, Luminosity

        and can be referred to either by these names or their SI unit names,
        e.g. length, metre, and m all refer to the same thing here.
        """
        if len(args):
            if isSequenceType(args[0]) and len(args[0]) == 7:
                # initialisation by list
                self._dims = args[0]
            elif isinstance(args[0], Dimension):
                # initialisation by another dimension object
                self._dims = args[0]._dims
        else:
            # initialisation by keywords
            self._dims = [0, 0, 0, 0, 0, 0, 0]
            for k in keywords.keys():
                # _di stores the index of the dimension with name 'k'
                self._dims[_di[k]] = keywords[k]
    #### METHODS ####
    def get_dimension(self, d):
        """Returns the list of dimension indices.

        See documentation for __init__.
        """
        return self._dims[_di[d]]

    def set_dimension(self, d, value):
        """Sets the list of dimension indices.

        See documentation for __init__.
        """
        self._dims[_di[d]] = value

    def is_dimensionless(self):
        """Tells you whether the object is dimensionless."""
        return sum([x == 0 for x in self._dims]) == 7
    #### REPRESENTATION ####
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        """String representation in basic SI units, or 1 for dimensionless."""
        s = ""
        for i in range(len(self._dims)):
            if self._dims[i]:
                s += _ilabel[i]
                if self._dims[i] != 1: s += "^" + str(self._dims[i])
                s += " "
        if not len(s): return "1"
        return s.strip()
    #### ARITHMETIC ####
    # Note that none of the dimension arithmetic objects do sanity checking
    # on their inputs, although most will throw an exception if you pass the
    # wrong sort of input
    def __mul__(self, value):
        return Dimension([x + y for x, y in izip(self._dims, value._dims)])

    def __div__(self, value):
        return Dimension([x - y for x, y in izip(self._dims, value._dims)])

    def __truediv__(self, value):
        return self.__div__(value)

    def __pow__(self, value):
        return Dimension([x * value for x in self._dims])

    def __imul__(self, value):
        self._dims = [x + y for x, y in izip(self._dims, value._dims)]
        return self

    def __idiv__(self, value):
        self._dims = [x - y for x, y in izip(self._dims, value._dims)]
        return self

    def __itruediv__(self, value):
        return self.__idiv__(value)

    def __ipow__(self, value):
        self._dims = [x * value for x in self._dims]
        return self
    #### COMPARISON ####
    def __eq__(self, value):
        #return sum([x==y for x,y in izip(self._dims,value._dims)])==7
        return sum([is_within_absolute_tolerance(x, y) for x, y in izip(self._dims, value._dims)]) == 7

    def __ne__(self, value):
        return not self.__eq__(value)
    #### MAKE DIMENSION PICKABLE ####
    def __getstate__(self):
        return self._dims

    def __setstate__(self, state):
        self._dims = state


class Scale(object):
    """Stores the scale factor for each SI dimension.

    Probably would only be very rarely used by a user, but might
    conceivably be useful in certain circumstances.

    Methods:

    -- Initialisation by list of keywords
    -- scale_factor(dim) gives the overall scaling for a value in
       dimension dim
    -- unit_representation(dim) gives a string representation of
       the unit defined by the Scale object applied to dimension
       dim
    """
    __slots__ = ["scale"]

    def __init__(self, *args, **keywords):
        """Initialise by list of scales or keywords, see Dimension documentation

        e.g. Scale(length="m", time="u") =
             Scale(["m","","u","","","",""])
        corresponds to measuring the unit of length at the milli scale
        and the unit of time at the u scale.
        """
        self.scale = [ "", "", "", "", "", "", "" ]
        for k in keywords:
            self.scale[_di[k]] = keywords[k]

    def scale_factor(self, dim):
        """Returns the scaling factor for dimension dim

        For example, if the scale factor of length is milli, and the
        dimensions of dim are length^2 then the scale factor will be
        0.001^2.
        """
        sf = 1
        for s, i in izip(self.scale, dim._dims):
            if i: sf *= _siprefixes[s] ** i
        return sf

    def unit_representation(self, dim):
        """Returns a representation of the dimension dim at this scale

        For example, if the scale factor of length is milli, and the
        dimensions of dim are length^2 then this will return mm^2.
        """
        s = ""
        for i in range(7):
            if dim._dims[i]:
                s += self.scale[i] + _ilabel[i]
                if dim._dims[i] != 1:
                    s += "^" + str(dim._dims[i])
                s += " "
        return s.strip()


class DimensionMismatchError(Exception):
    """Exception class for attempted operations with inconsistent dimensions

    For example, ``3*mvolt + 2*amp`` raises this exception. The purpose of this
    class is to help catch errors based on incorrect units. The exception will
    print a representation of the dimensions of the two inconsistent objects
    that were operated on. If you want to check for inconsistent units in your
    code, do something like::

        try:
            ...
            your code here
            ...
        except DimensionMismatchError, inst:
            ...
            cleanup code here, e.g.
            print "Found dimension mismatch, details:", inst
            ...
    """
    def __init__(self, description, *dims):
        """Raise as DimensionMismatchError(desc,dim1,dim2,...)

        desc -- a description of the type of operation being performed, e.g.
                Addition, Multiplication, etc.
        dim -- the dimensions of the objects involved in the operation, any
               number of them is possible
        """
        self._dims = dims
        self.desc = description

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = self.desc + ", dimensions were "
        for d in self._dims:
            s += "(" + str(d) + ") "
        return s

def is_scalar_type(obj):
    """Tells you if the object is a 1d number type

    This function is mostly used internally by the module for
    argument type checking. A scalar type can be considered
    a dimensionless quantity (see the documentation for
    Quantity for more information).
    """
    return isNumberType(obj) and not isSequenceType(obj)

def get_dimensions(obj):
    """Returns the dimensions of any object that has them.

    Slightly more general than obj.get_dimensions() because it will return
    a new dimensionless Dimension() object if the object is of number type
    but not a Quantity (e.g. a float or int).
    """
    if isNumberType(obj) and not isinstance(obj, Quantity): return Dimension()
    return obj.get_dimensions()


def is_dimensionless(obj):
    """Tests if a scalar value is dimensionless or not, returns a ``bool``.

    Note that the syntax may change in later releases of Brian, with tighter
    integration of scalar and array valued quantities.
    """
    return get_dimensions(obj) == Dimension()


def have_same_dimensions(obj1, obj2):
    """Tests if two scalar values have the same dimensions, returns a ``bool``.

    Note that the syntax may change in later releases of Brian, with tighter
    integration of scalar and array valued quantities.
    """
    return get_dimensions(obj1) == get_dimensions(obj2)

def display_in_unit(x, u):
    """String representation of the object x in unit u.
    """
    if not have_same_dimensions(x, u):
        raise DimensionMismatchError("Non-matching unit for function display_in_unit", get_dimensions(x), get_dimensions(u))
    s = str(float(x / u)) + " "
    if not is_dimensionless(u):
        if isinstance(u, Unit):
            s += str(u)
        else:
            s += str(u.dim)
    return s.strip()

def quantity_with_dimensions(floatval, dims):
    return Quantity.with_dimensions(floatval, dims)


class Quantity(numpy.float64):
    """A number with an associated physical dimension.

    In most cases, it is not necessary to create a :class:`Quantity` object
    by hand, instead use the constant unit names ``second``, ``kilogram``,
    etc. The details of how :class:`Quantity` objects work is subject to
    change in future releases of Brian, as we plan to reimplement it
    in a more efficient manner, more tightly integrated with numpy. The
    following can be safely used:

    * :class:`Quantity`, this name will not change, and the usage
      ``isinstance(x,Quantity)`` should be safe.
    * The standard unit objects, ``second``, ``kilogram``, etc.
      documented in the main documentation will not be subject
      to change (as they are based on SI standardisation).
    * Scalar arithmetic will work with future implementations.
    """

    # This documentation is subject to change.
    """
    This is the main user class for the units module, although
    in most cases it is not necessary to initialise a new
    quantity by hand (see construction below for details).

    The Quantity class defines arithmetic operations which
    check for consistency of dimensions and raise the
    DimensionMismatchError exception if they are inconsistent.

    The class also defines default and other representations
    of a number for printing purposes.

    Typical usage:

    I = 3 * amp # I is a Quantity object
    R = 2 * ohm # same for R
    print I*R # displays "6 V"
    print (I*R).in_unit(mvolt) # displays "6000 mV"
    print (I*R)/mvolt # displays "6000"
    x = I + R # raises DimensionMismatchError

    See the documentation on the Unit class for more details
    about the available unit names like mvolt, etc.

    Casting rules:

    The three rules that define the casting operations for
    Quantity object are:

    1. Quantity op Quantity = Quantity
        - Performs dimension checking if appropriate
    2. Scalar op Quantity = Quantity
        - Assumes that the scalar is dimensionless
    3. other op Quantity = other
        - The Quantity object is downcast to a float

    Scalar types are 1 dimensional number types, including float, int, etc.
    but not array.

    The Quantity class is a derived class of float, so many other operations
    will also downcast to float. For example, sin(x) where x is a quantity
    will return sin(float(x)) without doing any dimension checking.

    Construction details:

    x = Quantity(value) returns a dimensionless object, you can then
        set the dimensions via x.set_dimensions(dim)

    x = Quantity.with_dimensions(value,dim) returns an object with
        floating point value value, and dimensions dim, see the
        documentation for Quantity.with_dimensions(...) for more.

    Static constructors:

    -- with_dimensions(dim)
    -- with_dimensions(keywords...)

    Methods:

    -- get_dimensions() return Dimension
    -- set_dimensions(dim)
    -- is_dimensionless() return boolean
    -- at_scale(scale) return string
    -- has_same_dimensions(other) return boolean
    -- in_unit(unit) return string
    -- in_best_unit() return string
    """
    __slots__ = ["dim"]
    #### CONSTRUCTION ####
    def __init__(self, value):
        """Initialises as dimensionless
        """
        super(Quantity, self).__init__()
        self.dim = Dimension()
    @staticmethod
    def with_dimensions(value, *args, **keywords):
        """Static method to create a Quantity object with dimensions

        Use as Quantity.with_dimensions(value,dim),
               Quantity.with_dimensions(value,dimlist) or
               Quantity.with_dimensions(value,keywords...)

        -- value is a float or other scalar type
        -- dim is a dimension object
        -- dimlist, keywords (see the Dimension constructor)

        e.g.

        x = Quantity.with_dimensions(2,Dimension(length=1))
        x = Quantity.with_dimensions(2,length=1)
        x = 2 * metre

        all define the same object.
        """
        x = Quantity(value)
        if len(args) and isinstance(args[0], Dimension):
            x.set_dimensions(args[0])
        else:
            x.set_dimensions(Dimension(*args, **keywords))
        return x
    #### METHODS ####
    def get_dimensions(self):
        """Returns the dimensions of this object
        """
        return self.dim

    def set_dimensions(self, dim):
        """Set the dimensions of this object
        """
        self.dim = dim

    def is_dimensionless(self):
        """Tells you whether this is a dimensionless object
        """
        return self.dim.is_dimensionless()

    def at_scale(self, scale):
        """Returns a string representation at given scale
        """
        return str(float(self) / scale.scale_factor(self.dim)) + " " + scale.unit_representation(self.dim)

    def has_same_dimensions(self, other):
        """Tells you if this object has the same dimensions as another.
        """
        return self.dim == get_dimensions(other)

    def in_unit(self, u):
        """String representation of the object in unit u.
        """
        if not self.has_same_dimensions(u):
            raise DimensionMismatchError("Non-matching unit for method in_unit", self.dim, u.dim)
        s = str(float(self / u)) + " "
        if not u.is_dimensionless():
            if isinstance(u, Unit):
                s += str(u)
            else:
                s += str(u.dim)
        return s.strip()

    def in_best_unit(self, *regs):
        """String representation of the object in the 'best unit'

        For more information, see the documentation for the UnitRegistry
        class. Essentially, this looks at the value of the quantity for
        all 'known' matching units (e.g. mvolt, namp, etc.) and returns
        the one with the most compact representation. Standard units are
        built in, but you can register new units for consideration.
        """
        u = _get_best_unit(self, *regs)
        return self.in_unit(u)
    #### METHODS (NUMERICAL) ####
    def sqrt(self):
        return self ** 0.5

    def log(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.log(float(self)), self.dim)
        raise DimensionMismatchError('log', self.dim)

    def exp(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.exp(float(self)), self.dim)
        raise DimensionMismatchError('exp', self.dim)

    def sin(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.sin(float(self)), self.dim)
        raise DimensionMismatchError('sin', self.dim)

    def cos(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.cos(float(self)), self.dim)
        raise DimensionMismatchError('cos', self.dim)

    def tan(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.tan(float(self)), self.dim)
        raise DimensionMismatchError('tan', self.dim)

    def asin(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.asin(float(self)), self.dim)
        raise DimensionMismatchError('asin', self.dim)

    def acos(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.acos(float(self)), self.dim)
        raise DimensionMismatchError('acos', self.dim)

    def atan(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.atan(float(self)), self.dim)
        raise DimensionMismatchError('atan', self.dim)
    arcsin = asin
    arccos = cos
    arctan = tan

    def sinh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.sinh(float(self)), self.dim)
        raise DimensionMismatchError('sinh', self.dim)

    def cosh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.cosh(float(self)), self.dim)
        raise DimensionMismatchError('cosh', self.dim)

    def tanh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(math.tanh(float(self)), self.dim)
        raise DimensionMismatchError('tanh', self.dim)

    def arcsinh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(numpy.arcsinh(float(self)), self.dim)
        raise DimensionMismatchError('sinh', self.dim)

    def arccosh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(numpy.arccosh(float(self)), self.dim)
        raise DimensionMismatchError('cosh', self.dim)

    def arctanh(self):
        if self.is_dimensionless():
            return Quantity.with_dimensions(numpy.arctanh(float(self)), self.dim)
        raise DimensionMismatchError('tanh', self.dim)
    #### REPRESENTATION ####
    def __repr__(self):
        #return self.in_best_unit()
        return self.__str__()

    def __str__(self):
        #s = super(Quantity,self).__str__()
        #if not self.is_dimensionless(): s += " " + str(self.dim)
        #return s
        return self.in_best_unit()
        #return str(float(self))+'*'+str(get_unit(self))
    #### ARITHMETIC ####
    # Arithmetic operations implement the following set of rules for
    # determining casting:
    # 1. Quantity op Quantity returns Quantity (and performs dimension checking if appropriate)
    # 2. Scalar op Quantity returns Quantity (and performs dimension checking assuming Scalar is dimensionless)
    # 3. other op Quantity returns other (Quantity is downcast to float)
    # Scalar types are those for which is_scalar_type() returns True, including float, int, long, complex but not array
    def __mul__(self, other):
        # This code, like all the other arithmetic code below, implements the casting rules
        # defined above.
        if isinstance(other, Quantity):
            return Quantity.with_dimensions(float(self) * float(other), self.dim * other.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(float(self) * other, self.dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__mul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        if isinstance(other, Quantity):
            return Quantity.with_dimensions(float(self) / float(other), self.dim / other.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(float(self) / other, self.dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__div__(other)

    def __truediv__(self, other):
        if isinstance(other, Quantity):
            return Quantity.with_dimensions(float(self) / float(other), self.dim / other.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(float(self) / other, self.dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__truediv__(other)

    def __rdiv__(self, other):
        if isinstance(other, Quantity):
            return Quantity.with_dimensions(float(other) / float(self), other.dim / self.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(other / float(self), [-x for x in self.dim._dims])
        else:
            return NotImplemented
            #return super(Quantity,self).__rdiv__(other)

    def __rtruediv__(self, other):
        if isinstance(other, Quantity):
            return Quantity.with_dimensions(float(other) / float(self), other.dim / self.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(other / float(self), [-x for x in self.dim._dims])
        else:
            return NotImplemented
            #return super(Quantity,self).__rtruediv__(other)

    def __mod__(self, other):
        if isinstance(other, Quantity) or is_scalar_type(other):
            dim = get_dimensions(other)
            if dim == self.dim:
                return Quantity.with_dimensions(float(self) % float(other), self.dim)
            else: raise DimensionMismatchError("Addition", self.dim, dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__add__(other)

    def __add__(self, other):
        if isinstance(other, Quantity) or is_scalar_type(other):
            dim = get_dimensions(other)
            if dim == self.dim:
                return Quantity.with_dimensions(float(self) + float(other), self.dim)
            else: raise DimensionMismatchError("Addition", self.dim, dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__add__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Quantity) or is_scalar_type(other):
            dim = get_dimensions(other)
            if dim == self.dim:
                return Quantity.with_dimensions(float(self) - float(other), self.dim)
            else: raise DimensionMismatchError("Subtraction", self.dim, dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__sub__(other)

    def __rsub__(self, other):
        if isinstance(other, Quantity) or is_scalar_type(other):
            dim = get_dimensions(other)
            if dim == self.dim:
                return Quantity.with_dimensions(float(other) - float(self), self.dim)
            else: raise DimensionMismatchError("Subtraction(R)", self.dim, dim)
        else:
            return NotImplemented
            #return super(Quantity,self).__rsub__(other)

    def __pow__(self, other):
        if isinstance(other, Quantity):
            if other.is_dimensionless():
                # WARNING: because dimension consistency is checked by exact comparison of dimensions,
                # this may lead to unexpected behaviour (e.g. (x**2)**0.5 may not have the same dimensions as x)
                return Quantity.with_dimensions(float(self) ** float(other), self.dim ** float(other))
            else: raise DimensionMismatchError("Power", self.dim, other.dim)
        elif is_scalar_type(other):
            return Quantity.with_dimensions(float(self) ** other, self.dim ** other)
        else:
            return NotImplemented
            #return super(Quantity,self).__pow__(other)

    def __rpow__(self, other):
        if self.is_dimensionless():
            if isinstance(other, Quantity):
                return Quantity.with_dimensions(float(other) ** float(self), other.dim ** float(self))
            elif is_scalar_type(other):
                return Quantity(other ** float(self))
            else:
                return NotImplemented
                #return super(Quantity,self).__pow__(other)
        else:
            raise DimensionMismatchError("Power(R)", self.dim)

    def __neg__(self):
        return Quantity.with_dimensions(-float(self), self.dim)

    def __pos__(self):
        return self

    def __abs__(self):
        return Quantity.with_dimensions(abs(float(self)), self.dim)
    #### COMPARISONS ####
    def __lt__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) < float(other)
            else: raise DimensionMismatchError("LessThan", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0.: return float(self) < other
            if numpy.isposinf(other): return True
            if numpy.isneginf(other): return False
            if self.is_dimensionless():
                return float(self) < other
            else: raise DimensionMismatchError("LessThan", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__lt__(other)

    def __le__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) <= float(other)
            else: raise DimensionMismatchError("LessThanOrEquals", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0.: return float(self) <= other
            if numpy.isposinf(other): return True
            if numpy.isneginf(other): return False
            if self.is_dimensionless():
                return float(self) <= other
            else: raise DimensionMismatchError("LessThanOrEquals", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__le__(other)

    def __gt__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) > float(other)
            else: raise DimensionMismatchError("GreaterThan", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0.: return float(self) > other
            if numpy.isneginf(other): return True
            if numpy.isposinf(other): return False
            if self.is_dimensionless():
                return float(self) > other
            else: raise DimensionMismatchError("GreaterThan", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__gt__(other)

    def __ge__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) >= float(other)
            else: raise DimensionMismatchError("GreaterThanOrEquals", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0.: return float(self) >= other
            if numpy.isneginf(other): return True
            if numpy.isposinf(other): return False
            if self.is_dimensionless():
                return float(self) >= other
            else: raise DimensionMismatchError("GreaterThanOrEquals", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__ge__(other)

    def __eq__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) == float(other)
            else: raise DimensionMismatchError("Equals", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0. or numpy.isinf(other): return float(self) == other
            if self.dim.is_dimensionless():
                return float(self) == other
            else: raise DimensionMismatchError("Equals", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__eq__(other)

    def __ne__(self, other):
        if isinstance(other, Quantity):
            if self.dim == other.dim:
                return float(self) != float(other)
            else: raise DimensionMismatchError("Equals", self.dim, other.dim)
        elif is_scalar_type(other):
            if other == 0 or other == 0. or numpy.isinf(other): return float(self) != other
            if self.dim.is_dimensionless():
                return float(self) != other
            else: raise DimensionMismatchError("NotEquals", self.dim, Dimension())
        else:
            return NotImplemented
            #return super(Quantity,self).__ne__(other)
    #### MAKE QUANTITY PICKABLE ####
    def __reduce__(self):
        return (quantity_with_dimensions, (float(self), self.dim))


class Unit(Quantity):
    '''
    A physical unit

    Normally, you do not need to worry about the implementation of
    units. They are derived from the :class:`Quantity` object with
    some additional information (name and string representation).
    You can define new units which will be used when generating
    string representations of quantities simply by doing an
    arithmetical operation with only units, for example::

        Nm = newton * metre

    Note that operations with units are slower than operations with
    :class:`Quantity` objects, so for efficiency if you do not need the
    extra information that a :class:`Unit` object carries around, write
    ``1*second`` in preference to ``second``.
    '''

    # original documentation
    """A physical unit

    Basically, a unit is just a quantity with given dimensions, e.g.
    mvolt = 0.001 with the dimensions of voltage. The units module
    defines a large number of standard units, and you can also define
    your own (see below).

    The unit class also keeps track of various things that were used
    to define it so as to generate a nice string representation of it.
    See Representation below.

    Typical usage:

    x = 3 * mvolt # returns a quantity
    print x.in_unit(uvolt) # returns 3000 uV

    Standard units:

    The units class has the following fundamental units:

    metre, kilogram, second, amp, kelvin, mole, candle

    And these additional basic units:

    radian, steradian, hertz, newton, pascal, joule, watt,
    coulomb, volt, farad, ohm, siemens, weber, tesla, henry,
    celsius, lumen, lux, becquerel, gray, sievert, katal

    And additionally, it includes all scaled versions of these
    units using the following prefixes

     Factor     Name    Prefix
     -----      ----    ------
     10^24      yotta   Y
     10^21      zetta   Z
     10^18      exa     E
     10^15      peta    P
     10^12      tera    T
     10^9       giga    G
     10^6       mega    M
     10^3       kilo    k
     10^2       hecto   h
     10^1       deka    da
     1
     10^-1      deci    d
     10^-2      centi   c
     10^-3      milli   m
     10^-6      micro   u (\mu in SI)
     10^-9      nano    n
     10^-12     pico    p
     10^-15     femto   f
     10^-18     atto    a
     10^-21     zepto   z
     10^-24     yocto   y

    So for example nohm, ytesla, etc. are all defined.

    Defining your own:

    It can be useful to define your own units for printing
    purposes. So for example, to define the newton metre, you
    write:

    Nm = newton * metre

    Writing:

    print (1*Nm).in_unit(Nm)

    will return "1 Nm" because the Unit class generates a new
    display name of "Nm" from the display names "N" and "m" for
    newtons and metres automatically (see Representation below).

    To register this unit for use in the automatic printing
    of the Quantity.in_best_unit() method, see the documentation
    for the UnitRegistry class.

    Construction:

    The best way to construct a new unit is to use standard units
    already defined and arithmetic operations, e.g. newton*metre.
    See the documentation for __init__ and the static methods create(...)
    and create_scaled_units(...) for more details.

    If you don't like the automatically generated display name for
    the unit, use the set_display_name(name) method.

    Representation:

    A new unit defined by multiplication, division or taking powers
    generates a name for the unit automatically, so that for
    example the name for pfarad/mmetre**2 is "pF/mm^2", etc. If you
    don't like the automatically generated name, use the
    set_display_name(name) method.
    """
    __slots__ = ["dim", "scale", "scalefactor", "dispname", "name", "iscompound"]
    #### CONSTRUCTION ####
    def __init__(self, value):
        """Initialises a dimensionless unit
        """
        super(Unit, self).__init__(value)
        self.dim = Dimension()
        self.scale = [ "", "", "", "", "", "", "" ]
        self.scalefactor = ""
        self.dispname = ""
        self.iscompound = False

    def __new__(typ, *args, **kw):
        obj = super(Unit, typ).__new__(typ, *args, **kw)
        global automatically_register_units
        if automatically_register_units:
            register_new_unit(obj)
        return obj
    @staticmethod
    def create(dim, name="", dispname="", scalefactor="", **keywords):
        """Creates a new named unit

        dim -- the dimensions of the unit
        name -- the full name of the unit, e.g. volt
        dispname -- the display name, e.g. V
        scalefactor -- scaling factor, e.g. m for mvolt
        keywords -- the scaling for each SI dimension, e.g. length="m", mass="-1", etc.
        """
        scale = [ "", "", "", "", "", "", "" ]
        for k in keywords:
            scale[_di[k]] = keywords[k]
        v = 1.0
        for s, i in izip(scale, dim._dims):
            if i: v *= _siprefixes[s] ** i
        u = Unit(v * _siprefixes[scalefactor])
        u.dim = dim
        u.scale = scale
        u.scalefactor = scalefactor + ""
        u.name = name + ""
        u.dispname = dispname + ""
        u.iscompound = False
        return u
    @staticmethod
    def create_scaled_unit(baseunit, scalefactor):
        """Create a scaled unit from a base unit

        baseunit -- e.g. volt, amp
        scalefactor -- e.g. "m" for mvolt, mamp
        """
        u = Unit(float(baseunit) * _siprefixes[scalefactor])
        u.dim = baseunit.dim
        u.scale = baseunit.scale
        u.scalefactor = scalefactor
        u.name = scalefactor + baseunit.name
        u.dispname = scalefactor + baseunit.dispname
        u.iscompound = False
        return u
    #### METHODS ####
    def set_name(self, name):
        """Sets the name for the unit
        """
        self.name = name

    def set_display_name(self, name):
        """Sets the display name for the unit
        """
        self.dispname = name
    #### REPRESENTATION ####
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        if self.dispname == "":
            s = self.scalefactor + " "
            for i in range(7):
                if self.dim._dims[i]:
                    s += self.scale[i] + _ilabel[i]
                    if self.dim._dims[i] != 1: s += "^" + str(self.dim._dims[i])
                    s += " "
            if not len(s): return "1"
            return s.strip()
        else:
            return self.dispname
    #### ARITHMETIC ####
    def __mul__(self, other):
        if isinstance(other, Unit):
            u = Unit(float(self) * float(other))
            u.name = self.name + other.name
            u.dispname = self.dispname + ' ' + other.dispname
            u.dim = self.dim * other.dim
            u.iscompound = True
            return u
        else:
            return super(Unit, self).__mul__(other)

    def __div__(self, other):
        if isinstance(other, Unit):
            u = Unit(float(self) / float(other))
            u.name = self.name + 'inv_' + other.name + '_endinv'
            if other.iscompound:
                u.dispname = '(' + self.dispname + ')'
            else:
                u.dispname = self.dispname
            u.dispname += '/'
            if other.iscompound:
                u.dispname += '(' + other.dispname + ')'
            else:
                u.dispname += other.dispname
            u.dim = self.dim / other.dim
            u.iscompound = True
            return u
        else:
            return super(Unit, self).__div__(other)

    def __pow__(self, other):
        if is_scalar_type(other):
            u = Unit(float(self) ** other)
            u.name = self.name + 'pow_' + str(other) + '_endpow'
            if self.iscompound:
                u.dispname = '(' + self.dispname + ')'
            else:
                u.dispname = self.dispname
            u.dispname += '^' + str(other)
            u.dim = self.dim ** other
            return u
        else:
            return super(Unit, self).__mul__(other)

automatically_register_units = False
#### FUNDAMENTAL UNITS
metre = Unit.create(Dimension(m=1), "metre", "m")
meter = Unit.create(Dimension(m=1), "meter", "m")
kilogram = Unit.create(Dimension(kg=1), "kilogram", "kg")
gram = Unit.create_scaled_unit(kilogram, "m")
gram.set_name('gram')
gram.set_display_name('g')
gramme = Unit.create_scaled_unit(kilogram, "m")
gramme.set_name('gramme')
gramme.set_display_name('g')
second = Unit.create(Dimension(s=1), "second", "s")
amp = Unit.create(Dimension(A=1), "amp", "A")
kelvin = Unit.create(Dimension(K=1), "kelvin", "K")
mole = Unit.create(Dimension(mol=1), "mole", "mol")
candle = Unit.create(Dimension(candle=1), "candle", "cd")
fundamental_units = [ metre, meter, gram, second, amp, kelvin, mole, candle ]

#### DERIVED UNITS, from http://physics.nist.gov/cuu/Units/units.html
derived_unit_table = \
        [\
        [ 'radian', 'rad', Dimension() ], \
        [ 'steradian', 'sr', Dimension() ], \
        [ 'hertz', 'Hz', Dimension(s= -1) ], \
        [ 'newton', 'N', Dimension(m=1, kg=1, s= -2) ], \
        [ 'pascal', 'Pa', Dimension(m= -1, kg=1, s= -2) ], \
        [ 'joule', 'J', Dimension(m=2, kg=1, s= -2) ], \
        [ 'watt', 'W', Dimension(m=2, kg=1, s= -3) ], \
        [ 'coulomb', 'C', Dimension(s=1, A=1) ], \
        [ 'volt', 'V', Dimension(m=2, kg=1, s= -3, A= -1) ], \
        [ 'farad', 'F', Dimension(m= -2, kg= -1, s=4, A=2) ], \
        [ 'ohm', 'ohm', Dimension(m=2, kg=1, s= -3, A= -2) ], \
        [ 'siemens', 'S', Dimension(m= -2, kg= -1, s=3, A=2) ], \
        [ 'weber', 'Wb', Dimension(m=2, kg=1, s= -2, A= -1) ], \
        [ 'tesla', 'T', Dimension(kg=1, s= -2, A= -1) ], \
        [ 'henry', 'H', Dimension(m=2, kg=1, s= -2, A= -2) ], \
        [ 'celsius', 'degC', Dimension(K=1) ], \
        [ 'lumen', 'lm', Dimension(cd=1) ], \
        [ 'lux', 'lx', Dimension(m= -2, cd=1) ], \
        [ 'becquerel', 'Bq', Dimension(s= -1) ], \
        [ 'gray', 'Gy', Dimension(m=2, s= -2) ], \
        [ 'sievert', 'Sv', Dimension(m=2, s= -2) ], \
        [ 'katal', 'kat', Dimension(s= -1, mol=1) ]\
        ]

# Pointless list only here so that static analysis in Eclipse works ok
# All the values here are overwritten by the code below
volt = Unit(1); mvolt = Unit(1); uvolt = Unit(1)
namp = Unit(1); mamp = Unit(1); uamp = Unit(1); pamp = Unit(1)
ohm = Unit(1); Mohm = Unit(1); kohm = Unit(1)
siemens = Unit(1); msiemens = Unit(1); usiemens = Unit(1)
hertz = Unit(1); khertz = Unit(1); Mhertz = Unit(1)
farad = Unit(1); mfarad = Unit(1); ufarad = Unit(1); nfarad = Unit(1)
msecond = Unit(1)

# Generate derived unit objects and make a table of base units from these and the fundamental ones
base_units = fundamental_units + [gramme, kilogram] # make a copy
for _du in derived_unit_table:
    _u = Unit.create(_du[2], _du[0], _du[1])
    exec _du[0] + "=_u"
    base_units.append(_u)

all_units = base_units + []

# Generate scaled units for all base units
scaled_units = []
for _bu in base_units:
    for _k in _siprefixes.keys():
        if len(_k):
            _u = Unit.create_scaled_unit(_bu, _k)
            exec _k + _bu.name + "=_u"
            all_units.append(_u)
            if not _k in ["da", "d", "c", "h"]:
                scaled_units.append(_u)

# Generate 2nd and 3rd powers for all scaled base units
powered_units = []
for bu in all_units + []:
    for i in [2, 3]:
        u = bu ** i
        u.name = bu.name + str(i)
        exec bu.name + str(i) + '=u'
        all_units.append(u)
        if not bu.scalefactor in ['da', 'd', 'c', 'h']:
            powered_units.append(u)

# Define additional units

# Current list from http://physics.nist.gov/cuu/Units/units.html, far from complete
additional_units = [ pascal * second, newton * metre, watt / metre ** 2, joule / kelvin, \
                   joule / (kilogram * kelvin), joule / kilogram, watt / (metre * kelvin), \
                   joule / metre ** 3, volt / metre ** 3, coulomb / metre ** 3, coulomb / metre ** 2, \
                   farad / metre, henry / metre, joule / mole, joule / (mole * kelvin), \
                   coulomb / kilogram, gray / second, katal / metre ** 3 ]

automatically_register_units = True


class UnitRegistry(object):
    """Stores known units for printing in best units

    All a user needs to do is to use the register_new_unit(u)
    function.

    Default registries:

    The units module defines three registries, the standard units,
    user units, and additional units. Finding best units is done
    by first checking standard, then user, then additional. New
    user units are added by using the register_new_unit(u) function.

    Standard units includes all the basic non-compound unit names
    built in to the module, including volt, amp, etc. Additional
    units defines some compound units like newton metre (Nm) etc.

    Methods:

    add(u) - add a new unit
    __getitem__(x) - get the best unit for quantity x
      e.g. UnitRegistry ur; ur[3*mvolt] returns mvolt
    """
    def __init__(self):
        self.objs = []

    def add(self, u):
        """Add a unit to the registry
        """
        self.objs.append(u)

    def __getitem__(self, x):
        """Returns the best unit for quantity x

        The algorithm is to consider the value:

        m=abs(x/u)

        for all matching units u. If there is a unit u with a value of
        m in [1,1000) then we select that unit. Otherwise, we select
        the first matching unit (which will typically be the unscaled
        version).
        """
        matching = filter(lambda o: have_same_dimensions(o, x), self.objs)
        if len(matching) == 0:
            raise KeyError("Unit not found in registry.")
        floatrep = filter(lambda o: 0.1 <= abs(float(x / o)) < 100, matching)
        if len(floatrep):
            return floatrep[0]
        else:
            return matching[0]

def register_new_unit(u):
    """Register a new unit for automatic displaying of quantities

    Example usage:

    2.0*farad/metre**2 = 2.0 m^-4 kg^-1 s^4 A^2
    register_new_unit(pfarad / mmetre**2)
    2.0*farad/metre**2 = 2000000.0 pF/mm^2
    """
    UserUnitRegister.add(u)

standard_unit_register = UnitRegistry()
additional_unit_register = UnitRegistry()
UserUnitRegister = UnitRegistry()
map(standard_unit_register.add, base_units + scaled_units + powered_units)
map(additional_unit_register.add, additional_units)

def all_registered_units(*regs):
    """Returns all registered units in the correct order
    """
    if not len(regs):
        regs = [ standard_unit_register, UserUnitRegister, additional_unit_register]
    for r in regs:
        for u in r.objs:
            yield u

def _get_best_unit(x, *regs):
    """Returns the best unit for quantity x

    Checks the registries regs, unless none are provided in which
    case it will check the standard, user and additional unit
    registers in turn.
    """
    if get_dimensions(x) == Dimension():
        return Quantity(1)
    if len(regs):
        for r in regs:
            try:
                return r[x]
            except KeyError:
                pass
        return Quantity.with_dimensions(1, x.dim)
    else:
        return _get_best_unit(x, standard_unit_register, UserUnitRegister, additional_unit_register)

def get_unit(x, *regs):
    '''
    Find the most appropriate consistent unit from the unit registries, or just return a Quantity with the same dimensions and value 1
    '''
    for u in all_registered_units(*regs):
        if is_equal(float(u), 1) and have_same_dimensions(u, x):
            return u
    return Quantity.with_dimensions(1, get_dimensions(x))

def get_unit_fast(x):
    '''
    Return a quantity with value 1 and the same dimensions
    '''
    return Quantity.with_dimensions(1, get_dimensions(x))

#### DECORATORS


def check_units(**au):
    """Decorator to check units of arguments passed to a function

    **Sample usage:** ::

        @check_units(I=amp,R=ohm,wibble=metre,result=volt)
        def getvoltage(I,R,**k):
            return I*R

    You don't have to check the units of every variable in the function, and
    you can define what the units should be for variables that aren't
    explicitly named in the definition of the function. For example, the code
    above checks that the variable wibble should be a length, so writing::

        getvoltage(1*amp,1*ohm,wibble=1)

    would fail, but::

        getvoltage(1*amp,1*ohm,wibble=1*metre)

    would pass.
    String arguments are not checked (e.g. ``getvoltage(wibble='hello')`` would pass).

    The special name ``result`` is for the return value of the function.

    An error in the input value raises a :exc:`DimensionMismatchError`, and an error
    in the return value raises an ``AssertionError`` (because it is a code
    problem rather than a value problem).

    **Notes**

    This decorator will destroy the signature of the original function, and
    replace it with the signature ``(*args, **kwds)``. Other decorators will
    do the same thing, and this decorator critically needs to know the signature
    of the function it is acting on, so it is important that it is the first
    decorator to act on a function. It cannot be used in combination with another
    decorator that also needs to know the signature of the function.
    """
    def do_check_units(f):
        @wraps(f)
        def new_f(*args, **kwds):
            newkeyset = kwds.copy()
            arg_names = f.func_code.co_varnames[0:f.func_code.co_argcount]
            for (n, v) in zip(arg_names, args[0:f.func_code.co_argcount]):
                newkeyset[n] = v
            for k in newkeyset.iterkeys():
                if (k in au.keys()) and not isinstance(newkeyset[k], str): # string variables are allowed to pass, the presumption is they name another variable
                    if not have_same_dimensions(newkeyset[k], au[k]):
                        raise DimensionMismatchError("Function " + f.__name__ + " variable " + k + " should have dimensions of " + str(au[k]), get_dimensions(newkeyset[k]))
            result = f(*args, **kwds)
            if "result" in au:
                assert have_same_dimensions(result, au["result"]), \
                    "Function " + f.__name__ + " should return a value with unit " + str(au["result"]) + " but has returned " + str(get_dimensions(result))
            return result
#        new_f.__name__ = f.__name__
#        new_f.__doc__ = f.__doc__
#        new_f.__dict__.update(f.__dict__)
        return new_f
    return do_check_units

# Note: do not normally call this, see note on importing of decorator module at the top of this module
if use_decorator:
    old_check_units = check_units
    def check_units(**au):
        return lambda f : decorator.new_wrapper(old_check_units(**au)(f), f)
    check_units.__doc__ = old_check_units.__doc__

def _check_nounits(**au):
    """Don't bother checking units decorator
    """
    def dont_check_units(f):
        return f
    return dont_check_units

def scalar_representation(x):
    if isinstance(x, Unit):
        return x.name
    u = get_unit(x)
    if isinstance(u, Unit):
        return '(' + repr(float(x)) + '*' + u.name + ')'
    if isinstance(x, Quantity):
        return '(Quantity.with_dimensions(' + repr(float(x)) + ',' + repr(x.dim._dims) + '))'
    return repr(x)

# Remove all units
if not bup.use_units:
    for _u in all_units:
        exec _u.name + "=float(_u)"
    check_units = _check_nounits
    def get_dimensions(obj):
        return Dimension()

    def is_dimensionless(obj):
        return True

    def have_same_dimensions(obj1, obj2):
        return True

    def get_unit(x, *regs):
        return 1.

    def scalar_representation(x):
        return '1.0'

# Add unit names to __all__
all_unit_names = [u.name for u in all_units]
__all__.extend(all_unit_names)

if __name__ == "__main__":

#    # the pattern 'pat' below is a regular expression for all the unit names
#    # you can use it as an exclusion pattern for the epydoc api docs
#    base_unit_ids = set([id(_) for _ in base_units])
#    l = [_ for _ in __all__ if id(locals()[_]) in base_unit_ids]
#    anybaseunit = '('+'|'.join(l)+')'
#    print anybaseunit
#    prefixes = [_ for _ in _siprefixes.keys() if _]
#    anyprefix = '('+'|'.join(prefixes)+')'
#    print anyprefix
#    pat = anyprefix+'?'+anybaseunit+'[23]?'
#    print pat
#    import re
#    for x in __all__:
#        if not len(re.findall(pat,x)):
#            print x

    from numpy import *

    # shorthand function used for example code below
    def pE(vname, str):
        uname = vname
        if vname == "": uname = "temp"
        exec(uname + "=" + str)
        if vname != "": print vname, "=",
        print str,
        if locals()[uname] != None:
            print '=', locals()[uname]
        else:
            print
        return locals()[uname]

    V = pE("V", "3 * volt")
    I = pE("I", "2 * amp")
    a = pE("a", "array([1,2,3])")
    print
    R = pE("R", "V/I")
    pE("", "I*R")
    print
    pE("", "a*V")
    pE("", "V*a")
    pE("", "a+V")
    print
    pE("", "1000*metre")
    pE("", "1000*mmetre")
    print
    pE("", "(2*volt).in_unit(mvolt)")
    pE("", "(2*volt)/mvolt")
    print "(2*volt).in_unit(amp) =",
    try:
        print (2 * volt).in_unit(amp)
    except DimensionMismatchError, i:
        print "DimensionMismatchError:", i
    print
    pE("", "have_same_dimensions(1*volt,1*amp*ohm)")
    pE("", "have_same_dimensions(1*volt*second,1*amp*ohm)")
    print
    pE("", "(ufarad/nmetre)**2")
    print
    pE("", "2.0*farad/metre**2")
    pE("", "register_new_unit(pfarad / mmetre**2)")
    pE("", "2.0*farad/metre**2")

    print
    print "Some decorator examples (see code):"
    print

    @check_units(I=amp, R=ohm, wibble=metre, result=volt)
    def getvoltage(I, R, *args, **k):
        return I * R

    try:
        print getvoltage(1 * amp, 2 * ohm, 20)
        print getvoltage(R=2 * ohm, I=1 * amp, wibble=7 * mmetre)
        print getvoltage(1 * amp, 2 * ohm * metre)
    except DimensionMismatchError, inst:
        print "DME:", inst

    print
    pE("", "get_unit(3*msecond)")

# To avoid problems with Sage classes and Units
for k, v in globals().items():
    if isinstance(v, Unit):
        exec '* = Quantity.with_dimensions(float(*), *.dim)'.replace('*', k)

###################################################
##### ADDITIONAL INFORMATION

#SI DIMENSIONS
#-------------
#Quantity               Unit      Symbol
#--------               ----      ------
#Length                 metre     m
#Mass                   kilogram  kg
#Time                   second    s
#Electric current       ampere    A
#Temperature            kelvin    K
#Quantity of substance  mole      mol
#Luminosity             candle    cd

# SI UNIT PREFIXES
# ----------------
# Factor     Name    Prefix
# -----      ----    ------
# 10^24      yotta   Y
# 10^21      zetta   Z
# 10^18      exa     E
# 10^15      peta    P
# 10^12      tera    T
# 10^9       giga    G
# 10^6       mega    M
# 10^3       kilo    k
# 10^2       hecto   h
# 10^1       deka    da
# 1
# 10^-1      deci    d
# 10^-2      centi   c
# 10^-3      milli   m
# 10^-6      micro   u (\mu in SI)
# 10^-9      nano    n
# 10^-12     pico    p
# 10^-15     femto   f
# 10^-18     atto    a
# 10^-21     zepto   z
# 10^-24     yocto   y
