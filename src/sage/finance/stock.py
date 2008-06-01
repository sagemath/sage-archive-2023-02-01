"""
Stock market price series

AUTHORS:
    -- William Stein
"""

import urllib
from sage.structure.all import Sequence

from time_series import TimeSeries


class Day:
    def __init__(self, date, open, high, low, close, volume):
        self.date = date
        self.open=float(open); self.high=float(high); self.low=float(low); self.close=float(close)
        self.volume=int(volume)
    def __repr__(self):
        return '%10s %4.2f %4.2f %4.2f %4.2f %10d'%(self.date, self.open, self.high,
                   self.low, self.close, self.volume)

class Stock:
    def __init__(self, symbol):
        self.symbol = symbol.upper()

    def __repr__(self):
        return "%s (%s)"%(self.symbol, self.market_value())

    def market_value(self):
        return float(self.yahoo()['price'])

    def yahoo(self):
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
        try:
            return self.__historical
        except AttributeError:
            pass
        symbol = self.symbol
        def get_data(exchange):
            return urllib.urlopen(
              'http://finance.google.com/finance/historical?q=%s:%s&output=csv&startdate=%s'%(
                 exchange, symbol.upper(),startdate)).read()
        R = get_data('NASDAQ')
        if "Bad Request" in R:
             R = get_data("NYSE")
        R = R.splitlines()
        headings = R[0].split(',')
        self.__historical = []
        try:
            for x in reversed(R[1:]):
                date, opn, high, low, close, volume = x.split(',')
                self.__historical.append(Day(date, opn,high,low,close,volume))
        except ValueError:
             pass
        self.__historical = Sequence(self.__historical,cr=True,universe=lambda x:x)
        return self.__historical

    def close(self, *args, **kwds):
        return TimeSeries([x.close for x in self.historical(*args, **kwds)])
