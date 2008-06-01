"""
Markov Switching Multfractal model

REFERENCE: How to Forecast Long-Run Volatility: Regime Switching and
the Estimation of Multifractal Processes, Calvet and Fisher, 2004.
"""
import math
import random

from time_series import TimeSeries

class MarkovSwitchingMultifractal:
    def __init__(self, kbar, m0, sigma, gamma_kbar, b):
        """
        INPUT:
            kbar   -- positive integer
            m0     -- float with 0 <= m0 <= 2
            sigma  -- positive float
            gamma_kbar -- float with 0 < gamma_kbar < 1
            b      -- float > 1

        EXAMPLES:
            sage: msm = finance.MarkovSwitchingMultifractal(8,1.4,0.5,0.95,3); msm
            Markov switching multifractal model with m0 = 1.4, sigma = 0.5, b = 3.0, and gamma_8 = 0.95
            sage: yen_usd = finance.MarkovSwitchingMultifractal(10,1.448,0.461,0.998,3.76)
            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11)
        """
        self.__m0 = float(m0)
        assert self.__m0 >= 0 and self.__m0 <= 2, "m0 must be between 0 and 2"
        self.__sigma = float(sigma)
        assert self.__sigma > 0, "sigma must be positive"
        self.__b = float(b)
        assert self.__b > 1, "b must be bigger than 1"
        self.__gamma_kbar = float(gamma_kbar)
        assert self.__gamma_kbar > 0 and self.__gamma_kbar < 1, \
               "gamma_kbar must be between 0 and 1"
        self.__kbar = int(kbar)
        assert self.__kbar > 0, "kbar must be positive"

    def __repr__(self):
        return "Markov switching multifractal model with m0 = %s, sigma = %s, b = %s, and gamma_%s = %s"%(self.m0(), self.sigma(), self.b(), self.kbar(), self.gamma_kbar())

    def m0(self):
        return self.__m0
    def sigma(self):
        return self.__sigma
    def b(self):
        return self.__b
    def gamma_kbar(self):
        return self.__gamma_kbar
    def kbar(self):
        return self.__kbar

    def gamma(self):
        """
        Return the vector of the kbar transitional probabilities.

        OUTPUT:
            gamma -- a tuple of self.kbar() floats

        EXAMPLES:
            sage: msm = finance.MarkovSwitchingMultifractal(1.4,0.5,3,0.95,8)
            sage: msm.gamma()
            (0.001368852970712986, 0.0041009402016725094, 0.012252436441829828, 0.036308782091905023, 0.10501923017634662, 0.28312883556311919, 0.63159685013597011, 0.95000000000000351)
        """
        try:
            return self.__gamma
        except AttributeError:
            pass

        b          = self.__b
        gamma_kbar = self.__gamma_kbar
        kbar       = self.__kbar

        # We compute gamma1 frm gamma_kbar by inverting the relation
        # that defines the gamma_k given on page 54 of Calvet-Fisher:
        gamma1 = 1 - math.exp(math.log(1-gamma_kbar)/(b**(kbar-1)))

        gamma  = tuple([1 - (1 - gamma1)**(b**k) for k in range(kbar)])
        self.__gamma = gamma
        return gamma

    def simulation(self, T):
        """
        Retun a time series of the T values taken by simulating the
        running of this Markov switching multifractal model for T time
        steps.

        INPUT:
            T -- a positive integer
        OUTPUT:
            list -- a list of floats, the returns of the
                    log-prices of a financial asset or
                    exchange rate.

        EXAMPLES:
            sage: cad_usd = finance.MarkovSwitchingMultifractal(10,1.278,0.262,0.644,2.11); cad_usd
            Markov switching multifractal model with m0 = 1.278, sigma = 0.262, b = 2.11, and gamma_10 = 0.644
            sage: v = cad_usd.simulation(100); v
            [0.0011, -0.0032, 0.0006, 0.0007, 0.0034 ... -0.0023, 0.0008, 0.0015, -0.0003, 0.0027]
            sage: v
            [0.0011, -0.0032, 0.0006, 0.0007, 0.0034 ... -0.0023, 0.0008, 0.0015, -0.0003, 0.0027]
            sage: v.sums()
            [0.0011, -0.0021, -0.0015, -0.0008, 0.0026 ... -0.0383, -0.0376, -0.0360, -0.0363, -0.0336]
            sage: v.sums().exp().plot()
        """
        # Two values of the distribution M.
        m0 = self.m0()
        vals = [m0, 2 - m0]
        def M():
            return random.choice(vals)

        # Initalize the Markov volatility state vector
        from sage.rings.all import RDF
        kbar = self.kbar()
        m = (RDF**kbar)([M() for _ in xrange(kbar)])

        sigma = self.sigma()/100.0

        # r is the time series of returns
        r = TimeSeries(T)

        # Generate T Gaussian random numbers with mean 0
        # and variance 1.
        import scipy.stats
        eps = scipy.stats.norm().rvs(T)

        # Generate uniform distribution between 0 and 1
        uniform = scipy.stats.uniform().rvs(kbar*T)

        # The gamma_k
        gamma = self.gamma()

        for t in range(T):
            r[t] = sigma * eps[t] * m.prod().sqrt()
            # Now update the volatility state vector
            j = t*kbar
            for k in range(kbar):
                if uniform[j+k] <= gamma[k]:
                    # Draw from the distribution
                    m[k] = M()

        return r





def ml_estimation(v, kbar, M):
    """
    Compute parameters that model the time series v,

    INPUT:
        v -- series of returns; e.g., sequence of
             differences of logs of price
        kbar -- positive integer; model parameter
        m -- finite list of the values that the multiplier
             M takes on.

    OUTPUT:
        m0, sigma, gamma_kbar, b
    """

