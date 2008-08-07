"""
Stock market price series

AUTHORS:
    -- William Stein, 2008
    -- Brett Nakayama, 2008
    -- Chris Swierczewski, 2008

TESTS:
    sage: ohlc = sage.finance.stock.OHLC('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
    sage: loads(dumps(ohlc)) == ohlc
    True
"""

import urllib
from sage.structure.all import Sequence
from datetime import date

class OHLC:
    def __init__(self, timestamp, open, high, low, close, volume):
        """
        Open, high, low, and close information for a stock. Also stores
        a timestamp for that data along with the volume.

        INPUT:
            timestamp -- string
            open, high, low, close -- float
            volume -- int

        EXAMPLES:
            sage: sage.finance.stock.OHLC('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
             18-Aug-04 100.01 104.06 95.96 100.34   22353092
        """
        self.timestamp = timestamp
        self.open=float(open); self.high=float(high); self.low=float(low); self.close=float(close)
        self.volume=int(volume)

    def __repr__(self):
        """
        Return string representation of stock OHLC data.

        EXAMPLES:
            sage: sage.finance.stock.OHLC('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092).__repr__()
            ' 18-Aug-04 100.01 104.06 95.96 100.34   22353092'
        """
        return '%10s %4.2f %4.2f %4.2f %4.2f %10d'%(self.timestamp, self.open, self.high,
                   self.low, self.close, self.volume)

    def __cmp__(self, other):
        """
        Compare self and other.

        EXAMPLES:
            sage: ohlc = sage.finance.stock.OHLC('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
            sage: ohlc2 = sage.finance.stock.OHLC('18-Aug-04', 101.01, 104.06, 95.96, 100.34, 22353092)
            sage: cmp(ohlc, ohlc2)
            -1
        """
        if not isinstance(other, OHLC):
            return cmp(type(self), type(other))
        return cmp((self.timestamp, self.open, self.high, self.low, self.close, self.volume),
                   (other.timestamp, other.open, other.high, other.low, other.close, other.volume))

class Stock:
    """
    Class for retrieval of stock market information.
    """
    def __init__(self, symbol, cid=''):
        """
        Create a Stock object. Optional initialization by cid: an identifier
        for each equity used by Google Finance.

        INPUT:
            symbol -- string, a ticker symbol (with or without market)
                      format: "MARKET:SYMBOL" or "SYMBOL", if you don't
                              supply the market, it is assumed to be
                              NYSE or NASDAQ.
                      eg. "goog" or "OTC:NTDOY"
            cid    -- Integer, a Google contract ID (optional)

        LIMITATIONS:
            Currently, the symbol and cid do not have to match.  When using
            google(), the cid will take precedence.

        EXAMPLES:
            sage: S = finance.Stock('ibm')
            sage: S        # optional -- requires internet, and random
            IBM (127.48)
        """
        self.symbol = symbol.upper()
        self.cid = cid

    def __repr__(self):
        """
        Return string representation of this stock.

        EXAMPLES:
            sage: finance.Stock('ibm').__repr__()     # optional -- requires internet, and random
            'IBM (127.47)'
        """
        return "%s (%s)"%(self.symbol, self.market_value())

    def market_value(self):
        """
        Return the current market value of this stock.

        OUTPUT:
            Python float

        EXAMPLES:
            sage: finance.Stock('goog').market_value()   # optional and random
            575.83000000000004
        """
        return float(self.yahoo()['price'])

    def yahoo(self):
        """
        Get Yahoo current price data for this stock.

        OUTPUT:
            dict

        EXAMPLES:
            sage: finance.Stock('GOOG').yahoo()          # random and optional (requires internet)
            {'stock_exchange': '"NasdaqNM"', 'market_cap': '181.1B', '200day_moving_avg': '564.569', '52_week_high': '747.24', 'price_earnings_growth_ratio': '1.04', 'price_sales_ratio': '10.16', 'price': '576.48', 'earnings_per_share': '14.463', '50day_moving_avg': '549.293', 'avg_daily_volume': '6292480', 'volume': '1613507', '52_week_low': '412.11', 'short_ratio': '1.00', 'price_earnings_ratio': '40.50', 'dividend_yield': 'N/A', 'dividend_per_share': '0.00', 'price_book_ratio': '7.55', 'ebitda': '6.513B', 'change': '-9.32', 'book_value': '77.576'}
        """
        url = 'http://finance.yahoo.com/d/quotes.csv?s=%s&f=%s' % (self.symbol, 'l1c1va2xj1b4j4dyekjm3m4rr5p5p6s7')
        values = urllib.urlopen(url).read().strip().strip('"').split(',')
        data = {}
        data['price'] = values[0]
        data['change'] = values[1]
        data['volume'] = values[2]
        data['avg_daily_volume'] = values[3]
        data['stock_exchange'] = values[4]
        data['market_cap'] = values[5]
        data['book_value'] = values[6]
        data['ebitda'] = values[7]
        data['dividend_per_share'] = values[8]
        data['dividend_yield'] = values[9]
        data['earnings_per_share'] = values[10]
        data['52_week_high'] = values[11]
        data['52_week_low'] = values[12]
        data['50day_moving_avg'] = values[13]
        data['200day_moving_avg'] = values[14]
        data['price_earnings_ratio'] = values[15]
        data['price_earnings_growth_ratio'] = values[16]
        data['price_sales_ratio'] = values[17]
        data['price_book_ratio'] = values[18]
        data['short_ratio'] = values[19]
        return data

    def google(self,startdate='Jan+1,+1900',enddate=date.today().strftime("%b+%d,+%Y"), histperiod='daily'):
        """
        Return an immutable sequence of historical price data
        for this stock, obtained from Google. OHLC data is stored
        internally as well. By default, returns the past year's daily
        OHLC data.

        Dates startdate and enddate should be formatted 'Mon+d,+yyyy' where
        'Mon' is a three character abbreviation of the month's name.

        NOTE:
            Google Finance returns the past year's financial data by default
            when startdate is set too low from the equity's date of going
            public.  By default, this function only looks at the NASDAQ and
            NYSE markets.  However, if you sepcified the market during
            initialization of the stock (i.e. "finance.Stock("OTC:NTDOY")"),
            Stock.google() will give correct results.

        INPUT:
            startdate -- string, (default: 'Jan+1,+1900')
            enddate   -- string, (default: current date )
            histperiod -- string, ('daily' or 'weekly')

        OUTPUT:
            Sequence

        EXAMPLES:
        We get the first five days of VMware's stock history:
            sage: finance.Stock('vmw').google()[:5]   # optional -- requires internet
            [
             15-Aug-07 52.11 59.87 51.50 57.71   10678500,
             16-Aug-07 60.99 61.49 52.71 56.99    6919500,
             17-Aug-07 59.00 59.00 54.45 55.55    3086100,
             20-Aug-07 56.05 57.50 55.61 57.33    2140900,
             21-Aug-07 57.25 66.59 56.50 65.99    7369700
            ]

            sage: finance.Stock('F').google('Jan+3,+1978', 'Jul+7,+2008')[:5] # optional -- requires internet
            [
              3-Jan-78 2.90 2.90 2.84 2.85    1074495,
              4-Jan-78 2.84 2.84 2.81 2.83    1648713,
              5-Jan-78 2.83 2.84 2.77 2.77    1988524,
              6-Jan-78 2.77 2.77 2.74 2.76    2019988,
              9-Jan-78 2.73 2.73 2.69 2.73    2600499
            ]

        Note that when startdate is too far prior to a stock's actual start
        date, Google Finance defaults to a year's worth of stock history
        leading up to the specified enddate.  For example, Apple's (AAPL) stock
        history only dates back to September 7, 1984

            sage: finance.Stock('AAPL').google('Sep+1,+1900', 'Jan+1,+2000')[0:5] # optional -- requires internet
            [
              4-Jan-99 10.53 10.56 10.00 10.31   34031600,
              5-Jan-99 10.48 10.98 10.38 10.83   50361200,
              6-Jan-99 11.03 11.03 10.25 10.44   48163200,
              7-Jan-99 10.56 11.27 10.53 11.25   51036400,
              8-Jan-99 11.64 11.72 11.00 11.25   24244000
            ]

        Here is an example where we create and get the history of a stock
        that is not in NASDAQ or NYSE

            sage: finance.Stock("OTC:NTDOY").google(startdate="Jan+1,+2007", enddate="Jan+1,+2008")[:5]
            [
              3-Jan-07 32.44 32.75 32.30 32.44     156300,
              4-Jan-07 31.70 32.40 31.20 31.70     222700,
              5-Jan-07 30.15 30.50 30.15 30.15      65700,
              8-Jan-07 30.10 30.50 30.00 30.10     130800,
              9-Jan-07 29.90 30.05 29.60 29.90     103400
            ]


        Here, we create a stock by cid, and get historical data.
        Note that when using historical, if a cid is specified,
        it will take precedence over the stock's symbol.  So, if
        the symbol and cid do not match, the history based on the
        contract id will be returned.

            sage: sage.finance.stock.Stock("AAPL", 22144).google(startdate='Jan+1,+1990')[:5] #optional -- requires internet
            [
              2-Jan-90 8.81 9.38 8.75 9.31    6542800,
              3-Jan-90 9.50 9.50 9.38 9.38    7428400,
              4-Jan-90 9.56 9.69 9.31 9.41    7911200,
              5-Jan-90 9.44 9.56 9.25 9.44    4404000,
              8-Jan-90 9.38 9.50 9.25 9.50    3627600
            ]

        """
        cid = self.cid
        symbol = self.symbol

        if self.cid=='':
            if ':' in symbol:
                R = self._get_data('', startdate, enddate, histperiod)
            else:
                R = self._get_data('NASDAQ:', startdate, enddate, histperiod)
                if "Bad Request" in R:
                     R = self._get_data("NYSE:", startdate, enddate, histperiod)
        else:
            R = self._get_data('', startdate, enddate, histperiod)
        if "Bad Request" in R:
            raise RuntimeError
        self.__historical = []
        self.__historical = self._load_from_csv(R)
        return self.__historical

    def open(self, *args, **kwds):
        r"""
        Return a TimeSeries containing historical opening prices for this
        stock. If no arguments are given, will return last acquired historical
        data. Otherwise, data will be gotten from Google Finance.

        INPUT:
            startdate  -- string, (default: 'Jan+1,+1900')
            enddate    -- string, (default: current date)
            histperiod -- string, ('daily' or 'weekly')

        OUTPUT:
            TimeSeries -- Close price data

        EXAMPLES:
        You can directly obtain Open data as so:
            sage: finance.Stock('vmw').open(startdate='Jan+1,+2008', enddate='Feb+1,+2008')                 # optional -- requires internet
            [83.0500, 85.4900, 84.9000, 82.0000, 81.2500 ... 82.0000, 58.0500, 54.4900, 55.6000, 56.9800]

        Or, you can initialize stock data first and then extract the Open
        data:
            sage: c = finance.Stock('vmw')
            sage: c.google(startdate='Feb+1,+2008', enddate='Mar+1,+2008')[:5]    # optional -- requires internet
            [
             31-Jan-08 55.60 57.35 55.52 56.67    2607800,
              1-Feb-08 56.98 58.14 55.06 57.85    2489400,
              4-Feb-08 58.00 60.47 56.91 58.05    1840300,
              5-Feb-08 57.60 59.30 57.17 59.30    1711700,
              6-Feb-08 60.32 62.00 59.50 61.52    2209700
            ]
            sage: c.open()    # optional -- requires internet
            [55.6000, 56.9800, 58.0000, 57.6000, 60.3200 ... 56.5500, 59.3000, 60.0000, 59.7900, 59.2600]

        Otherwise, \code{self.google()} will be called with the default
        arguements returning a year's worth of data:
            sage: finance.Stock('vmw').open()   # optional and random -- requires internet and depends on day
            [52.1100, 60.9900, 59.0000, 56.0500, 57.2500 ... 83.0500, 85.4900, 84.9000, 82.0000, 81.2500]

        """

        from time_series import TimeSeries

        if len(args) != 0:
            return TimeSeries([x.open for x in self.google(*args, **kwds)])

        try:
            return TimeSeries([x.open for x in self.__historical])
        except AttributeError:
            pass

        return TimeSeries([x.open for x in self.google(*args, **kwds)])

    def close(self, *args, **kwds):
        r"""
        Return the time series of all historical closing prices for this stock.
        If no arguments are given, will return last aquired historical data.
        Otherwise, data will be gotten from Google Finance.

        INPUT:
            startdate  -- string, (default: 'Jan+1,+1900')
            enddate    -- string, (default: current date)
            histperiod -- string, ('daily' or 'weekly')

        OUTPUT:
            TimeSeries -- Close price data

        EXAMPLES:
        You can directly obtain close data as so:
            sage: finance.Stock('vmw').close(startdate='Jan+1,+2008', enddate='Feb+1,+2008')                 # optional -- requires internet
            [84.9900, 84.6000, 83.9500, 80.4900, 72.9900 ... 83.0000, 54.8700, 56.4200, 56.6700, 57.8500]

        Or, you can initialize stock data first and then extract the Close
        data:
            sage: c = finance.Stock('vmw')
            sage: c.google(startdate='Feb+1,+2008', enddate='Mar+1,+2008')[:5]    # optional -- requires internet
            [
             31-Jan-08 55.60 57.35 55.52 56.67    2607800,
              1-Feb-08 56.98 58.14 55.06 57.85    2489400,
              4-Feb-08 58.00 60.47 56.91 58.05    1840300,
              5-Feb-08 57.60 59.30 57.17 59.30    1711700,
              6-Feb-08 60.32 62.00 59.50 61.52    2209700
            ]
            sage: c.close()    # optional -- requires internet
            [56.6700, 57.8500, 58.0500, 59.3000, 61.5200 ... 58.2900, 60.1800, 59.8600, 59.9500, 58.6700]



        Otherwise, \code{self.google()} will be called with the default
        arguements returning a year's worh of data:
            sage: finance.Stock('vmw').close()   # optional and random -- requires internet and depends on day
            [57.7100, 56.9900, 55.5500, 57.3300, 65.9900 ... 84.9900, 84.6000, 83.9500, 80.4900, 72.9900]
        """

        from time_series import TimeSeries

        if len(args) != 0:
            return TimeSeries([x.close for x in self.google(*args, **kwds)])

        try:
            return TimeSeries([x.close for x in self.__historical])
        except AttributeError:
            pass

        return TimeSeries([x.close for x in self.google(*args, **kwds)])

    def load_from_file(self, file):
        """
        Load historical data from a local csv formatted data file. Note
        that no symbol data is included in Google Finance's csv data.
        The csv file must be formatted in the following way, just as
        on Google Finance:

        Timestamp,Open,High,Low,Close,Volume

        INPUT:
            file -- local file with Google Finance formatted OHLC data

        OUTPUTS:
            Sequence -- OHLC data

        EXAMPLES:
        Suppose you have a file in your home directoy containing Apple stock
        OHLC data, such as that from Google Finance, called AAPL-minutely.csv.
        One can load this information into a Stock object like so. Note that
        the path must be explicit:
            sage: finance.Stock('aapl').load_from_file(SAGE_ROOT + '/examples/finance/AAPL-minutely.csv')[:5]
            [
            1212408060 188.00 188.00 188.00 188.00        687,
            1212408000 188.00 188.11 188.00 188.00       2877,
            1212407700 188.00 188.00 188.00 188.00       1000,
            1212407640 187.75 188.00 187.75 188.00       2000,
            1212405780 187.80 187.80 187.80 187.80        100
            ]


        Note that since the source file doesn't contain information on which
        equity the information comes from, the symbol designated at
        initialization of Stock need not match the source of the data. For
        example, we can initialize a Stock object with the symbol 'goog',
        but load data from 'aapl' stock prices:
            sage: finance.Stock('goog').load_from_file(SAGE_ROOT + '/examples/finance/AAPL-minutely.csv')[:5]
            [
            1212408060 188.00 188.00 188.00 188.00        687,
            1212408000 188.00 188.11 188.00 188.00       2877,
            1212407700 188.00 188.00 188.00 188.00       1000,
            1212407640 187.75 188.00 187.75 188.00       2000,
            1212405780 187.80 187.80 187.80 187.80        100
            ]

        This tests a file that doesn't exist:
            sage: finance.Stock("AAPL").load_from_file("I am not a file")
            Bad path or file name

        """
        try:
            file_obj = open(file, 'r')
            R = file_obj.read();
            self.__historical = self._load_from_csv(R)
            file_obj.close()
            return self.__historical
        except IOError, msg:
            print "Bad path or file name"

    def _load_from_csv(self, R):
        """
        EXAMPLES:
        This indirectly tests _load_from_csv():
            sage: finance.Stock('aapl').load_from_file(SAGE_ROOT + "/examples/finance/AAPL-minutely.csv")
            [
            1212408060 188.00 188.00 188.00 188.00        687,
            1212408000 188.00 188.11 188.00 188.00       2877,
            1212407700 188.00 188.00 188.00 188.00       1000,
            1212407640 187.75 188.00 187.75 188.00       2000,
            1212405780 187.80 187.80 187.80 187.80        100
            ]


        """
        R = R.splitlines()
        headings = R[0].split(',')
        hist_data = []
        for x in reversed(R[1:]):
            try:
                timestamp, opn, high, low, close, volume = x.split(',')
                ohlc = OHLC(timestamp, opn,high,low,close,volume)
                hist_data.append(ohlc)
            except ValueError:
                pass
        hist_data = Sequence(hist_data,cr=True,universe=lambda x:x, immutable=True)
        return hist_data

    def _get_data(self, exchange='', startdate='Jan+1,+1900', enddate=date.today().strftime("%b+%d,+%Y"), histperiod='daily'):
        """
        This function is used internally.

        EXAMPLES:
        This indirectly tests the use of get_data.
            sage: finance.Stock('aapl').google(startdate='Jan+1,+1990')[:2]    # optional -- requires internet
            [
              2-Jan-90 8.81 9.38 8.75 9.31    6542800,
              3-Jan-90 9.50 9.50 9.38 9.38    7428400
            ]
        """
        symbol = self.symbol
        cid = self.cid
        if cid == '':
            url = 'http://finance.google.com/finance/historical?q=%s%s&startdate=%s&enddate=%s&histperiod=%s&output=csv'%(exchange, symbol.upper(), startdate, enddate, histperiod)
        else:
            url = 'http://finance.google.com/finance/historical?cid=%s&startdate=%s&enddate=%s&histperiod=%s&output=csv'%(cid, startdate, enddate, histperiod)
        return urllib.urlopen(url).read()
