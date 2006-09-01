#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from math import *
from sage.misc.misc import srange


def tasteful_ticks(minval, maxval):
    # Initialize the domain flags:
    onlyneg, domneg, onlypos, dompos = False, False, False, False

    # Is plot domain: *only or dominantly* and  *negative or positive*?
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
	    d1 = eval(s[n])

	#the stepsize will be 1 or greater:
	else:
	    sl = [s for s in str(int(abs(minval)))]
	    p = len(sl) #p is power for the order: o0=10^(p-1),o1=10^(p-2)
	    #d0 is first digit, d1 is second
	    d0 = eval(sl[0])
	    d1 = 0
	    #if p>1 we have steps of 1,10,100...
	    if p > 1:
                d1 = eval(sl[1])

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
	    d1 = eval(s[n])
	#the stepsize will be 1 or greater:
	else:
	    sl = [s for s in str(int(abs(maxval)))]
	    p = len(sl)
	    d0 = eval(sl[0])
	    d1 = 0
	    if p > 1:
                d1 = eval(sl[1])

    #choose a step size depending on either
    #1st or 2nd digits, d0 and d1, of given maxval or minval
    o0 = 10**(p-1)
    o1 = 10**(p-2)
    #fundamental step sizes: [1,2,2.5,5,10]
    if d0 == 1:
	if d1 > 5 and p!=1:
	    funda = 2.5
	    step = funda*o1
	elif d1 < 5 and p==1:
	    funda = 2.5
	    step = funda*o1
        else:
	    funda = 2
	    step = funda*o1
    elif d0 == 2 or d0 == 3:
	funda = 5
        step = funda*o1
    elif d0 == 4 or d0 == 5 or d0 == 6:
	funda = 1
        step = funda*o0
    else:
	funda = 2
        step = funda*o0

    #the 'fundamental' range
    #srange(step, d0*o0 + d1*o1 + step, step)
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

    if (abs((maxval - minval)/float(max(abs(maxval),abs(minval)))) < 0.2):
        tslmajor, oppaxis, step = tasteless_ticks(minval, maxval, 10)
    else:
        tslmajor, oppaxis, step = tasteful_ticks(minval, maxval)

    tslminor = srange(minval+0.2*2*step, maxval+0.2*2*step+0.2*step, 0.2*step)  #-0.2*2*step
    tslminor.sort()
    z = tslmajor[0]
    eps = z-tslminor[0]
    for w in tslminor: #[1:]
        eps2 = z-w
        if abs(eps2) > abs(eps):
            break
    tslminor = [w + eps for w in tslminor]
    tslminor = list(set(tslminor).union([-x for x in tslminor]) )
    tslminor = [x for x in tslminor if x >= minval and x <= maxval]
    return oppaxis, step, tslminor, tslmajor





