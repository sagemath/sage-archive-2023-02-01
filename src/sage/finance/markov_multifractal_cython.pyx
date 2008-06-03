"""
Markov Switching Multfractal model

Cython code
"""

from sage.misc.randstate cimport randstate, current_randstate

cdef extern from "math.h":
    double sqrt(double)

from time_series cimport TimeSeries

def simulation(Py_ssize_t n, Py_ssize_t k,
               double m0, double sigma, double b,
               int kbar, gamma):
    cdef double m1 = 2 - m0
    cdef Py_ssize_t i, j, a, c
    cdef TimeSeries t, eps
    cdef TimeSeries markov_state_vector = TimeSeries(kbar)
    cdef TimeSeries gamma_vals = TimeSeries(gamma)
    cdef randstate rstate = current_randstate()

    sigma = sigma / 100.0  # model's sigma is a percent

    # output list of simulations
    S = []

    for i from 0 <= i < k:
        # Initalize the model
        for j from 0 <= j < kbar:
            markov_state_vector._values[j] = m0 if (rstate.c_random() & 1) else m1    # n & 1 means "is odd"
        t = TimeSeries(n)

        # Generate n normally distributed random numbers with mean 0
        # and variance 1.
        eps = TimeSeries(n)
        eps.randomize('normal')

        for a from 0 <= a < n:
            # Compute next step in the simulation
            t._values[a] = sigma * eps._values[a] * sqrt(markov_state_vector.prod())

            # Now update the volatility state vector
            j = a * kbar
            for c from 0 <= c < kbar:
                if rstate.c_rand_double() <= gamma_vals._values[c]:
                    markov_state_vector._values[k] = m0 if (rstate.c_random() & 1) else m1

        S.append(t)

    return S







