"""
Markov Switching Multfractal model

Cython code
"""

cdef extern from "stdlib.h":
    long random()

from time_series cimport TimeSeries

def simulation(Py_ssize_t n, Py_ssize_t k,
               double m0, double sigma, double b,
               double gamma_kbar, int kbar):
    cdef double m1 = 2 - m0
    cdef Py_ssize_t i, j
    cdef TimeSeries t
    cdef TimeSeries markov_state_vector = TimeSeries(kbar)

    sigma = sigma / 100.0  # model's sigma is a percent

    # output list of simulations
    S = []

    for i from 0 <= i < k:
        # Initalize the model
        for j from 0 <= j < kbar:
            markov_state_vector._values[j] = m0 if random()%2 else m1
        t = TimeSeries(n)






