"""
Gosper iterator

A class which serves as a stateful iterable for computing the terms of the continued fraction of `(a*x+b)/(c*x+d)`,
where `a, b, c, d` are integers, and `x` is a continued fraction.

EXAMPLES::

    sage: from sage.rings.continued_fraction_gosper import gosper_iterator
    sage: x = continued_fraction(pi)
    sage: it = iter(gosper_iterator(3,2,3,1,x))
    sage: Word(it, length='infinite')
    word: 1,10,2,2,1,4,1,1,1,97,4,1,2,1,2,45,6,4,9,1,27,2,6,1,4,2,3,1,3,1,15,2,1,1,2,1,1,2,32,1,...
"""
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR


class gosper_iterator(object):

    def __init__(self, a, b, c, d, x):
        """
        Construct the class.

        INPUT:

        - ``a, b, c, d`` -- Integer coefficients of the transformation.
        - ``x`` -- An instance of a continued fraction.

        OUTPUT:

        - The instance of gosper_iterator class.

        TESTS::

            sage: a = Integer(randint(-10,10)); b = Integer(randint(-10,10));
            sage: c = Integer(randint(-10,10)); d = Integer(randint(-10,10));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: x = continued_fraction(([1,2],[3,4])); i = iter(gosper_iterator(a,b,c,d,x))
            sage: l = list(i)
            sage: preperiod_length = i.output_preperiod_length
            sage: preperiod = l[:preperiod_length]
            sage: period = l[preperiod_length:]
            sage: continued_fraction((preperiod, period), x.value()) == continued_fraction((a*x.value()+b)/(c*x.value()+d))
            True
        """
        from sage.rings.continued_fraction import ContinuedFraction_periodic
        self.a = a
        self.b = b
        self.c = c
        self.d = d

        self.x = iter(x)

        self.states = set()
        self.states_to_currently_emitted = dict()

        self.currently_emitted = 0
        self.currently_read = 0

        # Rational or quadratic case
        if isinstance(x, ContinuedFraction_periodic):
            self.input_preperiod_length = x.preperiod_length()
            self.input_period_length = x.period_length()
        # Infinite case
        else:
            self.input_preperiod_length = +Infinity
            self.input_period_length = 0

        self.output_preperiod_length = 0

    def __iter__(self):
        """
        Return the iterable instance of the class.

        Is called upon `iter(gosper_iterator(a,b,c,d,x))`.

        TESTS::

            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: ig = iter(gosper_iterator(a,b,c,d,continued_fraction(pi))); icf = iter(continued_fraction((a*pi+b)/(c*pi+d)));
            sage: lg = [next(ig) for _ in range(10)]; lcf = [next(icf) for _ in range(10)];
            sage: lg == lcf
            True
        """
        return self

    def __next__(self):
        """
        Return the next term of the transformation.

        TESTS::

            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: ig = iter(gosper_iterator(a,b,c,d,continued_fraction(pi))); icf = iter(continued_fraction((a*pi+b)/(c*pi+d)));
            sage: for i in range(10):
            ....:     assert next(ig) == next(icf)
        """
        limit = 100
        while True:
            if self.currently_read >= self.input_preperiod_length:
                current_state = (
                    ('a', self.a),
                    ('b', self.b),
                    ('c', self.c),
                    ('d', self.d),
                    ('index', (self.currently_read - self.input_preperiod_length) % self.input_period_length)
                )
                # for state in self.states:
                #     if self.compare_dicts(state, current_state, ['currently_emitted']):
                #         self.output_preperiod_length = state['currently_emitted']
                #         raise StopIteration
                if current_state in self.states:
                    self.output_preperiod_length = self.states_to_currently_emitted[current_state]
                    raise StopIteration
                self.states.add(current_state)
                self.states_to_currently_emitted[current_state] = self.currently_emitted
                if len(self.states) > 100:
                    print("ERROR: Stopping iteration, danger of memory overflow.")
                    raise StopIteration

            if (self.c == 0 and self.d == 0):
                raise StopIteration

            ub = self.bound(self.a, self.c)
            lb = self.bound(self.a + self.b, self.c + self.d)
            s = -self.bound(self.c, self.d)

            if ub == lb and s < 1:
                self.emit(ub)
                return Integer(ub)
            else:
                self.ingest()

            limit -= 1
            if limit < 1:
                print("ERROR: Next loop iteration ran too many times.")
                raise StopIteration

    def emit(self, q):
        """
        Change the state of the iterator, emitting the term `q`.

        TESTS::

            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gi = gosper_iterator(a,b,c,d,continued_fraction(pi))
            sage: for i in range(10):
            ....:     gi.emit(i)
            sage: gi.currently_emitted
            10
        """
        self.currently_emitted += 1
        # This is being computed for the case when no states are being saved (still reading preperiod).
        if self.currently_read <= self.input_preperiod_length:
            self.output_preperiod_length = self.currently_emitted
        a = self.a
        b = self.b
        self.a = self.c
        self.b = self.d
        self.c = a - q * self.c
        self.d = b - q * self.d

    def ingest(self):
        """
        Change the state of the iterator, ingesting another term from the input continued fraction.

        TESTS::

            sage: a = Integer(randint(-100,100)); b = Integer(randint(-100,100));
            sage: c = Integer(randint(-100,100)); d = Integer(randint(-100,100));
            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gi = gosper_iterator(a,b,c,d,continued_fraction(pi))
            sage: for i in range(10):
            ....:     gi.ingest()
            sage: gi.currently_read
            10
        """
        try:
            p = next(self.x)
            self.currently_read += 1
            a = self.a
            c = self.c
            self.a = a * p + self.b
            self.b = a
            self.c = c * p + self.d
            self.d = c
        except StopIteration:
            self.b = self.a
            self.d = self.c

    @staticmethod
    def bound(n, d):
        """
        Helper function for division. Return infinity if denominator is zero.

        TESTS::

            sage: from sage.rings.continued_fraction_gosper import gosper_iterator
            sage: gosper_iterator.bound(1,0)
            +Infinity
        """
        if d == 0:
            return +Infinity
        else:
            return (n / d).floor()
