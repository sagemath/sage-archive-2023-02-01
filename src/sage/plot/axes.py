#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from math import floor, log
from sage.misc.misc import srange


def tasteful_ticks(minval, maxval):
    minval, maxval = float(minval), float(maxval)
    # Initialize the domain flags:
    onlyneg, domneg, onlypos, dompos = False, False, False, False

    # Is the range: *only or dominantly* and  *negative or positive*?
    if abs(minval) > abs(maxval):
        if maxval < 0:
            onlyneg = True
        else:
            domneg  = True

        # Is the stepsize going to be < 1?
        if abs(minval) < 1:
            n = 0
            s = str(minval).split('.')[1]
            for c in s:
                n+=1
                if c != '0':
                    break
            p = -(n-1)
            d0 = eval(s[n-1])
            #string may only be of length 1
            try:
                d1 = eval(s[n])
            except IndexError:
                d1 = 0

        #the stepsize will be 1 or greater:
        else:
            if abs(minval) >= 10:
                sl = [s for s in str(int(abs(minval)))]
                d0 = eval(sl[0])
                d1 = eval(sl[1])
                p = len(sl)
            else:
                sl = str(abs(minval)).split('.')
                d0 = eval(sl[0])
                d1 = eval(sl[1][0])
                p = 1

    else: #this means: abs(minval) < abs(maxval)
        if minval > 0:
             onlypos = True
        else:
             dompos = True
        #is the stepsize going to be < 1?
        if abs(maxval) < 1:
            n = 0
            s = str(maxval).split('.')[1]
            for c in s:
                n+=1
                if c != '0':
                    break
            p = -(n-1)
            d0 = eval(s[n-1])
            try:
                sn = s[n]
            except IndexError:
                sn = "0"
            d1 = eval(sn)
        #the stepsize will be 1 or greater:
        else:
            if maxval >= 10:
                sl = [s for s in str(int(abs(maxval)))]
                d0 = eval(sl[0])
                d1 = eval(sl[1])
                p = len(sl)
            else:
                sl = str(abs(maxval)).split('.')
                d0 = eval(sl[0])
                d1 = eval(sl[1][0])
                p = 1

    #choose a step size depending on either
    #1st or 2nd digits, d0 and d1, of given maxval or minval
    o0 = 10**(p-1)
    o1 = 10**(p-2)
    #fundamental step sizes: [1,2,2.5,5,10]
    # 'fundamental' means that we split the x and y
    # ranges into widths from the above step sizes
    if d0 == 1:
        if d1 > 5 and p!=1:
            funda = 2.5
            step = funda*o1
        elif d1 < 5 and p==1:
            funda = 2.5
            step = funda*o1
        elif d1 < 5 and p>1:
            funda = 2.5
            step = funda*o1
        else:
            funda = 5
            step = funda*o1
    elif d0 == 2 or (d0 == 3 and p > 1):
        funda = 5
        step = funda*o1
    elif d0 in [3, 4, 5, 6]:
        funda = 1
        step = funda*o0
    else:
        funda = 2
        step = funda*o0

    #the 'fundamental' range
    fundrange = srange(0, d0*o0 + d1*o1 + step, step)

    def in_range(v):
        return [x for x in v if x >= minval and x <= maxval]

    #now major axis range choices
    oppaxis = 0

    if onlyneg:
        tslmajor = in_range([-x for x in fundrange]) #tick step list (major ticks)
        tslmajor.sort()
        oppaxis = tslmajor[-1] #position the other axis at this value

    elif domneg:
        neg = [-x for x in fundrange]
        pos = [x for x in fundrange if x <= abs(maxval)]
        tslmajor = in_range(neg + pos)
        tslmajor.sort()

    elif dompos:
        pos = [x for x in fundrange]
        neg = [-x for x in fundrange if x <= abs(minval)]
        tslmajor = in_range(neg + pos)
        tslmajor.sort()

    else: #onlypos
        tslmajor = in_range(fundrange)
        tslmajor.sort()
        oppaxis = tslmajor[0]
    return tslmajor, oppaxis, step


def trunc(x, digits_before_the_decimal):
    s = '%f'%float(x)
    i = s.find('.')
    t = s[:i - digits_before_the_decimal]
    if digits_before_the_decimal > 0:
        t += '0'* digits_before_the_decimal
    return float(eval(t))


def tasteless_ticks(minval, maxval, num_pieces):
    minval0 = minval
    maxval0 = maxval
    rnd   = int(floor(log(maxval - minval)/log(10)))
    if rnd < 0:
        rnd -= 1

    step  = (maxval - minval)/float(num_pieces)
    minval = trunc(minval, rnd)
    maxval = trunc(maxval+step, rnd)

    step  = (maxval - minval)/float(num_pieces)
    tslmajor = srange(minval, minval+(num_pieces+1)*step, step)

    tslmajor = [x for x in tslmajor if x >= minval0 and x <= maxval0]

    oppaxis = 0
    if maxval <= 0:  # only negative
        oppaxis = tslmajor[-1]
    elif minval >= 0:
        oppaxis = tslmajor[0]

    return tslmajor, oppaxis, step

def find_axes(minval, maxval):
    """
    Try to find axis tick positions that are well spaced
    """
    if minval >= maxval:
        raise ValueError, "maxval >= minval is required"

    # If there is a small differences between max and min values
    # compared to the size of the largest of (abs(maxval), abs(minval))
    # the function 'tasteless_ticks' is used, which in common usage is rare.
    if (abs((maxval - minval)/float(max(abs(maxval),abs(minval)))) < 0.2):
        tslmajor, oppaxis, step = tasteless_ticks(minval, maxval, 10)
    else:
        tslmajor, oppaxis, step = tasteful_ticks(minval, maxval)

    min = tslmajor[0] - step
    tslminor = srange(min, maxval + 0.2*step, 0.2*step)
    tslminor = [x for x in tslminor if x >= minval and x <= maxval]
    return oppaxis, step, tslminor, tslmajor


