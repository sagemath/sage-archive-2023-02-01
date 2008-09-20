from stock import Stock

from markov_multifractal import MarkovSwitchingMultifractal

from time_series import TimeSeries, autoregressive_fit

from fractal import (stationary_gaussian_simulation,
                     fractional_gaussian_noise_simulation,
                     fractional_brownian_motion_simulation,
                     multifractal_cascade_random_walk_simulation)

## def TimeSeries(values):
##     """
##     Initialize new time series.

##     INPUT:
##         values -- integer (number of values) or an iterable of floats

##     EXAMPLES:
##     This implicity calls init.
##         sage: finance.TimeSeries([pi, 3, 18.2])
##         [3.1416, 3.0000, 18.2000]
##     """
##     # A dirty trick to avoid importing time_series (hence real_double_vector)
##     # every time sage starts up, but to make sure that the finance.TimeSeries
##     # function is VERY fast.  The first time this function is called it is
##     # immediately replaced by the fast compiled version in the time_series
##     # module.
##     from time_series import TimeSeries
##     globals()['TimeSeries'] = TimeSeries
##     return TimeSeries(values)
