class lazy_attribute(object):
    r"""
    A lazy attribute for an object is like a usual attribute, except
    that, instead of being computed when the object is constructed
    (i.e. in __init__), it is computed on the fly the first time it
    is accessed.

    For constant values attached to an object, lazy attributes provide
    a shorter syntax and automatic caching (unlike methods), while
    playing well with inheritance (like methods): a subclass can
    easily override a given attribute; you don't need to call the
    super class constructor, etc.

    Technically, a lazy_attribute is a non-data descriptor (see
    Invoking Descriptors in the Python reference manual).

    EXAMPLES:
    We create a class whose instances have a lazy attribute x
      sage: class A(object):
      ...       def __init__(self):
      ...           self.a=2 # just to have some data to calculate from
      ...
      ...       @lazy_attribute
      ...       def x(self):
      ...           print "calculating x in A"
      ...           return self.a + 1
      ...

    For an instance a of A, a.x is calculated the first time it is accessed,
    and then stored as a usual attribute:
      sage: a = A()
      sage: a.x
      calculating x in A
      3
      sage: a.x
      3

    We redo the same example, but opening the hood to see what happens to
    the internal dictionary of the object:
      sage: a = A()
      sage: a.__dict__
      {'a': 2}
      sage: a.x
      calculating x in A
      3
      sage: a.__dict__
      {'a': 2, 'x': 3}
      sage: a.x
      3
      sage: timeit('a.x') # random
      625 loops, best of 3: 89.6 ns per loop

    This shows that, after the first calculation, the attribute x
    becomes a usual attribute; in particular, there is no time penalty
    to access it.

    A lazy attribute may be set as usual, even before its first access,
    in which case the lazy calculation is completely ignored:

      sage: a = A()
      sage: a.x = 4
      sage: a.x
      4
      sage: a.__dict__
      {'a': 2, 'x': 4}


    Conditional definitions.

    The function calculating the attribute may return NotImplemented
    to declare that, after all, it is not able to do it. In that case,
    the attribute lookup proceeds in the super class hierarchy:

      sage: class B(A):
      ...       @lazy_attribute
      ...       def x(self):
      ...           if hasattr(self, "y"):
      ...               print "calculating x from y in B"
      ...               return self.y
      ...           else:
      ...               print "y not there; B does not define x"
      ...               return NotImplemented
      ...
      sage: b = B()
      sage: b.x
      y not there; B does not define x
      calculating x in A
      3
      sage: b = B()
      sage: b.y = 1
      sage: b.x
      calculating x from y in B
      1

    Attribute existence testing

    Testing for the existence of an attribute with hasattr currently
    always triggers its full calculation, which may not be desirable
    when the calculation is expensive:

      sage: a = A()
      sage: hasattr(a, "x")
      calculating x in A
      True

    It would be great if we could take over the control somehow, if at
    all possible without a special implementation of hasattr, so as to
    allow for something like:

      sage: class A (object):
      ...       @lazy_attribute
      ...       def x(self, existence_only=False):
      ...           if existence_only:
      ...               print "testing for x existence"
      ...               return True
      ...           else:
      ...               print "calculating x in A"
      ...               return 3
      ...
      sage: a = A()
      sage: hasattr(a, "x") # todo: not implemented
      testing for x existence
      sage: a.x
      calculating x in A
      3
      sage: a.x
      3

    Here is a full featured example, with both conditional definition
    and existence testing:

      sage: class B(A):
      ...       @lazy_attribute
      ...       def x(self, existence_only=False):
      ...           if hasattr(self, "y"):
      ...               if existence_only:
      ...                   print "testing for x existence in B"
      ...                   return True
      ...               else:
      ...                   print "calculating x from y in B"
      ...                   return self.y
      ...           else:
      ...               print "y not there; B does not define x"
      ...               return NotImplemented
      ...
      sage: b = B()
      sage: hasattr(b, "x") # todo: not implemented
      y not there; B does not define x
      testing for x existence
      True
      sage: b.x
      y not there; B does not define x
      calculating x in A
      3
      sage: b = B()
      sage: b.y = 1
      sage: hasattr(b, "x") # todo: not implemented
      testing for x existence in B
      True
      sage: b.x
      calculating x from y in B
      1

    TESTS:

    Old style and new style classes play a bit differently with
    @property and attribute setting:

      sage: class A:
      ...       @property
      ...       def x(self):
      ...           print "calculating x"
      ...           return 3
      ...
      sage: a = A()
      sage: a.x = 4
      sage: a.__dict__
      {'x': 4}
      sage: a.x
      4
      sage: a.__dict__['x']=5
      sage: a.x
      5

      sage: class A (object):
      ...       @property
      ...       def x(self):
      ...           print "calculating x"
      ...           return 3
      ...
      sage: a = A()
      sage: a.x = 4
      Traceback (most recent call last):
      ...
      AttributeError: can't set attribute
      sage: a.__dict__
      {}
      sage: a.x
      calculating x
      3
      sage: a.__dict__['x']=5
      sage: a.x
      calculating x
      3

    In particular, lazy_attributes need to be implemented as non-data
    descriptors for new style classes, so as to leave access to
    setattr. We now check that this implementation also works for old
    style classes (conditional definition does not work yet).

      sage: class A:
      ...       def __init__(self):
      ...           self.a=2 # just to have some data to calculate from
      ...
      ...       @lazy_attribute
      ...       def x(self):
      ...           print "calculating x"
      ...           return self.a + 1
      ...
      sage: a = A()
      sage: a.__dict__
      {'a': 2}
      sage: a.x
      calculating x
      3
      sage: a.__dict__
      {'a': 2, 'x': 3}
      sage: a.x
      3
      sage: timeit('a.x') # random
      625 loops, best of 3: 115 ns per loop

      sage: a = A()
      sage: a.x = 4
      sage: a.x
      4
      sage: a.__dict__
      {'a': 2, 'x': 4}

      sage: class B(A):
      ...       @lazy_attribute
      ...       def x(self):
      ...           if hasattr(self, "y"):
      ...               print "calculating x from y in B"
      ...               return self.y
      ...           else:
      ...               print "y not there; B does not define x"
      ...               return NotImplemented
      ...
      sage: b = B()
      sage: b.x                         # todo: not implemented
      y not there; B does not define x
      calculating x in A
      3
      sage: b = B()
      sage: b.y = 1
      sage: b.x
      calculating x from y in B
      1
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, a, A):
        result = self.f(a)
        if result is NotImplemented:
            return getattr(super(A, a),self.f.func_name)
        setattr(a, self.f.func_name, result)
        return result
