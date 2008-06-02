"""
Stock market price series

AUTHORS:
    -- William Stein, 2008

TESTS:
    sage: day = sage.finance.stock.Day('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
    sage: loads(dumps(day)) == day
    True
"""

import urllib
from sage.structure.all import Sequence

from time_series import TimeSeries


class Day:
    def __init__(self, date, open, high, low, close, volume):
        """
        Summary stats for a day that a stock traded.

        INPUT:
            date -- string
            open, high, low, close -- float
            volume -- int

        EXAMPLES:
            sage: sage.finance.stock.Day('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
             18-Aug-04 100.01 104.06 95.96 100.34   22353092
        """
        self.date = date
        self.open=float(open); self.high=float(high); self.low=float(low); self.close=float(close)
        self.volume=int(volume)

    def __repr__(self):
        """
        Return string representation of stock day.

        EXAMPLES:
            sage: sage.finance.stock.Day('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092).__repr__()
            ' 18-Aug-04 100.01 104.06 95.96 100.34   22353092'
        """
        return '%10s %4.2f %4.2f %4.2f %4.2f %10d'%(self.date, self.open, self.high,
                   self.low, self.close, self.volume)

    def __cmp__(self, other):
        """
        Compare self and other.

        EXAMPLES:
            sage: day = sage.finance.stock.Day('18-Aug-04', 100.01, 104.06, 95.96, 100.34, 22353092)
            sage: day2 = sage.finance.stock.Day('18-Aug-04', 101.01, 104.06, 95.96, 100.34, 22353092)
            sage: cmp(day, day2)
            -1
        """
        if not isinstance(other, Day):
            return cmp(type(self), type(other))
        return cmp((self.date, self.open, self.high, self.low, self.close, self.volume),
                   (other.date, other.open, other.high, other.low, other.close, other.volume))

class Stock:
    def __init__(self, symbol):
        """
        Create a Stock object.

        INPUT:
            symbol -- string, a ticker symbol

        EXAMPLES:
            sage: S = finance.Stock('ibm')
            sage: S        # optional -- requires internet, and random
            IBM (127.48)
        """
        self.symbol = symbol.upper()

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

    def historical(self,startdate='Jan+1,+1990'):
        """
        Return an immutable sequence of historical price data
        for this stock, obtained from Google.

        INPUT:
            startdate -- string, (default: 'Jan+1,+1990')

        OUTPUT:
            Sequence

        EXAMPLES:
        We get the first five days of VMware's stock history:
            sage: finance.Stock('vmw').historical()[:5]   # optional -- requires internet
            [
             15-Aug-07 52.11 59.87 51.50 57.71   10678500,
             16-Aug-07 60.99 61.49 52.71 56.99    6919500,
             17-Aug-07 59.00 59.00 54.45 55.55    3086100,
             20-Aug-07 56.05 57.50 55.61 57.33    2140900,
             21-Aug-07 57.25 66.59 56.50 65.99    7369700
            ]
        """
        try:
            return self.__historical
        except AttributeError:
            pass
        symbol = self.symbol
        def get_data(exchange=''):
            """
            This function is used internally.
            """
            url = 'http://finance.google.com/finance/historical?q=%s%s&output=csv&startdate=%s'%(
                 exchange, symbol.upper(),startdate)
            return urllib.urlopen(url).read()
        if ':' in symbol:
            R = get_data()
        else:
            R = get_data('NASDAQ:')
            if "Bad Request" in R:
                 R = get_data("NYSE:")
        if "Bad Request" in R:
            raise RuntimeError
        R = R.splitlines()
        headings = R[0].split(',')
        self.__historical = []
        try:
            for x in reversed(R[1:]):
                date, opn, high, low, close, volume = x.split(',')
                self.__historical.append(Day(date, opn,high,low,close,volume))
        except ValueError:
             pass
        self.__historical = Sequence(self.__historical,cr=True,universe=lambda x:x, immutable=True)
        return self.__historical

    def close(self, *args, **kwds):
        """
        Return the time series of all historical closing prices for this stock.

        EXAMPLES:
            sage: finance.Stock('vmw').close()                 # optional -- requires internet
            [57.7100, 56.9900, 55.5500, 57.3300, 65.9900 ...
        """
        return TimeSeries([x.close for x in self.historical(*args, **kwds)])
