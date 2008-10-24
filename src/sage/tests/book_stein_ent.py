# -*- coding: utf-8 -*-

"""
This file contains all the example code from the published book
'Elementary Number Theory: Primes, Congruences, and Secretes' by
William Stein, Springer-Verlag, 2009.
"""


from sage.all_cmdline import *;
import sage.plot.plot; sage.plot.plot.DOCTEST_MODE=True

def warning_function(f):
    import warnings

    def doctest_showwarning(message, category, filename, lineno, file=f):
        try:
            file.write(warnings.formatwarning(message, category, 'doctest', lineno))
        except IOError:
            pass # the file (probably stdout) is invalid
    return doctest_showwarning

def change_warning_output(file):
    import warnings
    warnings.showwarning = warning_function(file)
def example_0():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> prime_range(Integer(10),Integer(50))###line 262:_sage_    : prime_range(10,50)
[11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
"""

def example_1():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> [n for n in range(Integer(10),Integer(30)) if not is_prime(n)]###line 267:_sage_    : [n for n in range(10,30) if not is_prime(n)]
[10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28]
"""

def example_2():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> gcd(Integer(97),Integer(100))###line 485:_sage_    : gcd(97,100)
1
>>> gcd(Integer(97) * Integer(10)**Integer(15), Integer(19)**Integer(20) * Integer(97)**Integer(2))###line 487:_sage_    : gcd(97 * 10^15, 19^20 * 97^2)
97
"""

def example_3():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> factor(Integer(1275))###line 606:_sage_    : factor(1275)
3 * 5^2 * 17
>>> factor(Integer(2007))###line 608:_sage_    : factor(2007)
3^2 * 223
>>> factor(Integer(31415926535898))###line 610:_sage_    : factor(31415926535898)
2 * 3 * 53 * 73 * 2531 * 534697
"""

def example_4():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> n = Integer(74037563479561712828046796097429573142593188889231289084936232638972765034028266276891996419625117843995894330502127585370118968098286733173273108930900552505116877063299072396380786710086096962537934650563796359)###line 721:_sage_    : n = 74037563479561712828046796097429573142593188889231289084936232638972765034028266276891996419625117843995894330502127585370118968098286733173273108930900552505116877063299072396380786710086096962537934650563796359
>>> len(n.str(Integer(2)))###line 722:_sage_    : len(n.str(2))
704
>>> len(n.str(Integer(10)))###line 724:_sage_    : len(n.str(10))
212
>>> n.is_prime()              # this is instant###line 726:_sage_    : n.is_prime()              # this is instant
False
"""

def example_5():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> p = Integer(2)**Integer(32582657) - Integer(1)###line 932:_sage_    : p = 2^32582657 - 1
>>> p.ndigits()###line 933:_sage_    : p.ndigits()
9808358

\noindent{}Next we convert $p$ to a decimal string and look at some of the digits.



9808358

'12457502601536945540'

'11752880154053967871'
"""

def example_6():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> prime_pi(Integer(6))###line 1081:_sage_    : prime_pi(6)
3
>>> prime_pi(Integer(100))###line 1083:_sage_    : prime_pi(100)
25
>>> prime_pi(Integer(3000000))###line 1085:_sage_    : prime_pi(3000000)
216816
"""

def example_7():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> plot(prime_pi, Integer(1),Integer(1000), rgbcolor=(Integer(0),Integer(0),Integer(1)))###line 1090:_sage_    : plot(prime_pi, 1,1000, rgbcolor=(0,0,1))
"""

def example_8():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> P = plot(Li, Integer(2),Integer(10000), rgbcolor='purple')###line 1220:_sage_    : P = plot(Li, 2,10000, rgbcolor='purple')
>>> Q = plot(prime_pi, Integer(2),Integer(10000), rgbcolor='black')###line 1221:_sage_    : Q = plot(prime_pi, 2,10000, rgbcolor='black')
>>> R = plot(sqrt(x)*log(x),Integer(2),Integer(10000),rgbcolor='red')###line 1222:_sage_    : R = plot(sqrt(x)*log(x),2,10000,rgbcolor='red')
>>> show(P+Q+R,xmin=Integer(0), figsize=[Integer(8),Integer(3)])###line 1223:_sage_    : show(P+Q+R,xmin=0, figsize=[8,3])
"""

def example_9():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = Integers(Integer(3))###line 1447:_sage_    : R = Integers(3)
>>> list(R)###line 1448:_sage_    : list(R)
[0, 1, 2]
"""

def example_10():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = Integers(Integer(10))###line 1669:_sage_    : R = Integers(10)
>>> a = R(Integer(3))                 # create an element of Z/10Z###line 1670:_sage_    : a = R(3)                 # create an element of Z/10Z
>>> a.multiplicative_order()###line 1671:_sage_    : a.multiplicative_order()
4

Notice that the powers of $a$ are periodic with period $4$, i.e.,
there are four powers and they repeat:

>>> [a**i for i in range(Integer(15))]###line 1678:_sage_    : [a^i for i in range(15)]
[1, 3, 9, 7, 1, 3, 9, 7, 1, 3, 9, 7, 1, 3, 9]
"""

def example_11():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> euler_phi(Integer(2007))###line 1712:_sage_    : euler_phi(2007)
1332
"""

def example_12():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> n = Integer(20)###line 1761:_sage_    : n = 20
>>> k = euler_phi(n); k###line 1762:_sage_    : k = euler_phi(n); k
8
>>> [Mod(x,n)**k for x in range(n) if gcd(x,n) == Integer(1)]###line 1764:_sage_    : [Mod(x,n)^k for x in range(n) if gcd(x,n) == 1]
[1, 1, 1, 1, 1, 1, 1, 1]
"""

def example_13():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> for n in range(Integer(1),Integer(10)):###line 1843:_sage_    : for n in range(1,10):
...    print n, factorial(n-Integer(1)) % n, -Integer(1) % n
1 0 0
2 1 1
3 2 2
4 2 3
5 4 4
6 0 5
7 6 6
8 0 7
9 0 8
"""

def example_14():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> CRT(Integer(2),Integer(3), Integer(3), Integer(5))###line 1957:_sage_    : CRT(2,3, 3, 5)
-7
"""

def example_15():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> CRT_list([Integer(2),Integer(3),Integer(2)], [Integer(3),Integer(5),Integer(7)])###line 1964:_sage_    : CRT_list([2,3,2], [3,5,7])
23
"""

def example_16():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> xgcd(Integer(5),Integer(7))###line 2147:_sage_    : xgcd(5,7)
(1, -4, 3)
>>> xgcd(Integer(130),Integer(61))###line 2149:_sage_    : xgcd(130,61)
(1, 23, -49)
"""

def example_17():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> a = Mod(Integer(17), Integer(61))###line 2214:_sage_    : a = Mod(17, 61)
>>> a**(-Integer(1))###line 2215:_sage_    : a^(-1)
18
"""

def example_18():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> Integer(100).str(Integer(2))###line 2266:_sage_    : 100.str(2)
'1100100'
"""

def example_19():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> Integer(0)*Integer(2)**Integer(0) + Integer(0)*Integer(2)**Integer(1) + Integer(1)*Integer(2)**Integer(2) + Integer(0)*Integer(2)**Integer(3) + Integer(0)*Integer(2)**Integer(4) + Integer(1)*Integer(2)**Integer(5) + Integer(1)*Integer(2)**Integer(6)###line 2271:_sage_    : 0*2^0 + 0*2^1 + 1*2^2 + 0*2^3 + 0*2^4 + 1*2^5 + 1*2^6
100
"""

def example_20():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> Mod(Integer(7),Integer(100))**Integer(91)###line 2344:_sage_    : Mod(7,100)^91
43
"""

def example_21():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> Integer(7)**Integer(91)###line 2350:_sage_    : 7^91
80153343160247310515380886994816022539378033762994852007501964604841680190743
"""

def example_22():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> n = Integer(95468093486093450983409583409850934850938459083)###line 2429:_sage_    : n = 95468093486093450983409583409850934850938459083
>>> Mod(Integer(2),n)**(n-Integer(1))###line 2430:_sage_    : Mod(2,n)^(n-1)
34173444139265553870830266378598407069248687241

Note that factoring $n$ actually takes much longer than the above
computation (which was essentially instant).

>>> factor(n)                # takes up to a few seconds.###line 2437:_sage_    : factor(n)                # takes up to a few seconds.
1610302526747 * 59285812386415488446397191791023889
"""

def example_23():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> n = Integer(95468093486093450983409583409850934850938459083)###line 2522:_sage_    : n = 95468093486093450983409583409850934850938459083
>>> is_prime(n)###line 2523:_sage_    : is_prime(n)
False
"""

def example_24():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> for p in primes(Integer(100)):###line 2529:_sage_    : for p in primes(100):
...   if is_prime(Integer(2)**p - Integer(1)):
...       print p, Integer(2)**p - Integer(1)
2 3
3 7
5 31
7 127
13 8191
17 131071
19 524287
31 2147483647
61 2305843009213693951
89 618970019642690137449562111
"""

def example_25():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def is_prime_lucas_lehmer(p):###line 2552:_sage_    : def is_prime_lucas_lehmer(p):
...    s = Mod(Integer(4), Integer(2)**p - Integer(1))
...    for i in range(Integer(3), p+Integer(1)):
...        s = s**Integer(2) - Integer(2)
...    return s == Integer(0)
>>> # Check primality of 2^9941 - 1###line 2557:_sage_    : # Check primality of 2^9941 - 1
>>> is_prime_lucas_lehmer(Integer(9941))###line 2558:_sage_    : is_prime_lucas_lehmer(9941)
True
>>> # Check primality of 2^next_prime(1000)-1###line 2560:_sage_    : # Check primality of 2^next_prime(1000)-1
>>> is_prime_lucas_lehmer(next_prime(Integer(1000)))###line 2561:_sage_    : is_prime_lucas_lehmer(next_prime(1000))
False
"""

def example_26():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> for p in primes(Integer(20)):###line 2608:_sage_    : for p in primes(20):
...    print p, primitive_root(p)
2 1
3 2
5 2
7 3
11 2
13 2
17 3
19 2
"""

def example_27():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = PolynomialRing(Integers(Integer(13)),names=('x',)); (x,) = R._first_ngens(Integer(1))###line 2658:_sage_    : R.<x> = PolynomialRing(Integers(13))
>>> f = x**Integer(15) + Integer(1)###line 2659:_sage_    : f = x^15 + 1
>>> f.roots()###line 2660:_sage_    : f.roots()
[(12, 1), (10, 1), (4, 1)]
>>> f(Integer(12))###line 2662:_sage_    : f(12)
0
"""

def example_28():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = PolynomialRing(Integers(Integer(13)),names=('x',)); (x,) = R._first_ngens(Integer(1))###line 2696:_sage_    : R.<x> = PolynomialRing(Integers(13))
>>> f = x**Integer(6) + Integer(1)###line 2697:_sage_    : f = x^6 + 1
>>> f.roots()###line 2698:_sage_    : f.roots()
[(11, 1), (8, 1), (7, 1), (6, 1), (5, 1), (2, 1)]
"""

def example_29():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> log(RealNumber('19683.0'))###line 3348:_sage_    : log(19683.0)
9.88751059801299
>>> log(RealNumber('3.0'))###line 3350:_sage_    : log(3.0)
1.09861228866811
>>> log(RealNumber('19683.0')) / log(RealNumber('3.0'))###line 3352:_sage_    : log(19683.0) / log(3.0)
9.00000000000000
"""

def example_30():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> plot(log, RealNumber('0.1'),Integer(10), rgbcolor=(Integer(0),Integer(0),Integer(1)))###line 3401:_sage_    : plot(log, 0.1,10, rgbcolor=(0,0,1))
"""

def example_31():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> p = Integer(53)###line 3406:_sage_    : p = 53
>>> R = Integers(p)###line 3407:_sage_    : R = Integers(p)
>>> a = R.multiplicative_generator()###line 3408:_sage_    : a = R.multiplicative_generator()
>>> v = sorted([(a**n, n) for n in range(p-Integer(1))])###line 3409:_sage_    : v = sorted([(a^n, n) for n in range(p-1)])
>>> G = plot(point(v,pointsize=Integer(50),rgbcolor=(Integer(0),Integer(0),Integer(1))))###line 3410:_sage_    : G = plot(point(v,pointsize=50,rgbcolor=(0,0,1)))
>>> H = plot(line(v,rgbcolor=(RealNumber('0.5'),RealNumber('0.5'),RealNumber('0.5'))))###line 3411:_sage_    : H = plot(line(v,rgbcolor=(0.5,0.5,0.5)))
>>> G + H###line 3412:_sage_    : G + H
"""

def example_32():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> q = Integer(93450983094850938450983409623)###line 3504:_sage_    : q = 93450983094850938450983409623
>>> q.is_prime()###line 3505:_sage_    : q.is_prime()
True
>>> is_prime((q-Integer(1))//Integer(2))###line 3507:_sage_    : is_prime((q-1)//2)
True
>>> g = Mod(-Integer(2), q)###line 3509:_sage_    : g = Mod(-2, q)
>>> g.multiplicative_order()###line 3510:_sage_    : g.multiplicative_order()
93450983094850938450983409622
>>> n = Integer(18319922375531859171613379181)###line 3512:_sage_    : n = 18319922375531859171613379181
>>> m = Integer(82335836243866695680141440300)###line 3513:_sage_    : m = 82335836243866695680141440300
>>> g**n###line 3514:_sage_    : g^n
45416776270485369791375944998
>>> g**m###line 3516:_sage_    : g^m
15048074151770884271824225393
>>> (g**n)**m###line 3518:_sage_    : (g^n)^m
85771409470770521212346739540
>>> (g**m)**n###line 3520:_sage_    : (g^m)^n
85771409470770521212346739540
"""

def example_33():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def rsa(bits):###line 3671:_sage_    : def rsa(bits):
...    # only prove correctness up to 1024 bits
...    proof = (bits <= Integer(1024))
...    p = next_prime(ZZ.random_element(Integer(2)**(bits//Integer(2) +Integer(1))),
...                proof=proof)
...    q = next_prime(ZZ.random_element(Integer(2)**(bits//Integer(2) +Integer(1))),
...                proof=proof)
...    n = p * q
...    phi_n = (p-Integer(1)) * (q-Integer(1))
...    while True:
...        e = ZZ.random_element(Integer(1),phi_n)
...        if gcd(e,phi_n) == Integer(1): break
...    d = lift(Mod(e,phi_n)**(-Integer(1)))
...    return e, d, n
...
>>> def encrypt(m,e,n):###line 3686:_sage_    : def encrypt(m,e,n):
...    return lift(Mod(m,n)**e)
...
>>> def decrypt(c,d,n):###line 3689:_sage_    : def decrypt(c,d,n):
...    return lift(Mod(c,n)**d)
...
>>> e,d,n = rsa(Integer(20))###line 3692:_sage_    : e,d,n = rsa(20)
>>> c = encrypt(Integer(123), e, n)###line 3693:_sage_    : c = encrypt(123, e, n)
>>> decrypt(c, d, n)###line 3694:_sage_    : decrypt(c, d, n)
123
"""

def example_34():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def encode(s):###line 3758:_sage_    : def encode(s):
...    s = str(s)      # make input a string
...    return sum(ord(s[i])*Integer(256)**i for i in range(len(s)))
>>> def decode(n):###line 3761:_sage_    : def decode(n):
...    n = Integer(n)  # make input an integer
...    v = []
...    while n != Integer(0):
...        v.append(chr(n % Integer(256)))
...        n //= Integer(256)    # this replaces n by floor(n/256).
...    return ''.join(v)
>>> m = encode('Run Nikita!'); m###line 3768:_sage_    : m = encode('Run Nikita!'); m
40354769014714649421968722
>>> decode(m)###line 3770:_sage_    : decode(m)
'Run Nikita!'
"""

def example_35():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def crack_rsa(n, phi_n):###line 3902:_sage_    : def crack_rsa(n, phi_n):
...    R = PolynomialRing(QQ,names=('x',)); (x,) = R._first_ngens(Integer(1))
...    f = x**Integer(2) - (n+Integer(1) -phi_n)*x + n
...    return [b for b, _ in f.roots()]
>>> crack_rsa(Integer(31615577110997599711), Integer(31615577098574867424))###line 3906:_sage_    : crack_rsa(31615577110997599711, 31615577098574867424)
[8850588049, 3572144239]
"""

def example_36():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def crack_when_pq_close(n):###line 3959:_sage_    : def crack_when_pq_close(n):
...    t = Integer(ceil(sqrt(n)))
...    while True:
...       k = t**Integer(2) - n
...       if k > Integer(0):
...          s = Integer(int(round(sqrt(t**Integer(2) - n))))
...          if s**Integer(2) + n == t**Integer(2):
...             return t+s, t-s
...
...       t += Integer(1)
...
>>> crack_when_pq_close(Integer(23360947609))###line 3970:_sage_    : crack_when_pq_close(23360947609)
(153649, 152041)


For example, you might think that choosing a random
prime, and the next prime after would be a good idea,
but instead it creates an easy-to-crack cryptosystem.

>>> p = next_prime(Integer(2)**Integer(128)); p###line 3979:_sage_    : p = next_prime(2^128); p
340282366920938463463374607431768211507
>>> q = next_prime(p)###line 3981:_sage_    : q = next_prime(p)
>>> crack_when_pq_close(p*q)###line 3982:_sage_    : crack_when_pq_close(p*q)
(340282366920938463463374607431768211537,
      340282366920938463463374607431768211507)
"""

def example_37():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def crack_given_decrypt(n, m):###line 4125:_sage_    : def crack_given_decrypt(n, m):
...    n = Integer(n); m = Integer(m);  # some type checking
...    # Step 1: divide out powers of 2
...    while True:
...        if is_odd(m): break
...        divide_out = True
...        for i in range(Integer(5)):
...           a = randrange(Integer(1),n)
...           if gcd(a,n) == Integer(1):
...              if Mod(a,n)**(m//Integer(2)) != Integer(1):
...                  divide_out = False
...                  break
...        if divide_out:
...            m = m//Integer(2)
...        else:
...            break
...    # Step 2: Compute GCD
...    while True:
...        a = randrange(Integer(1),n)
...        g = gcd(lift(Mod(a, n)**(m//Integer(2))) - Integer(1), n)
...        if g != Integer(1) and g != n:
...            return g
...

\noindent{}We show how to verify Example~\ref{ex:crackprob} using \sage.

>>> n=Integer(32295194023343); e=Integer(29468811804857); d=Integer(11127763319273)###line 4152:_sage_    : n=32295194023343; e=29468811804857; d=11127763319273
>>> crack_given_decrypt(n, e*d - Integer(1))###line 4153:_sage_    : crack_given_decrypt(n, e*d - 1)
737531
>>> factor(n)###line 4155:_sage_    : factor(n)
737531 * 43788253


\noindent{}We try a much larger example.

>>> e = Integer(22601762315966221465875845336488389513)###line 4162:_sage_    : e = 22601762315966221465875845336488389513
>>> d = Integer(31940292321834506197902778067109010093)###line 4163:_sage_    : d = 31940292321834506197902778067109010093
>>> n = Integer(268494924039590992469444675130990465673)###line 4164:_sage_    : n = 268494924039590992469444675130990465673
>>> p = crack_given_decrypt(n, e*d - Integer(1))###line 4165:_sage_    : p = crack_given_decrypt(n, e*d - 1)
>>> print "ignore this";  p   # random output (could be other prime divisor)###line 4166:_sage_    : p   # random output (could be other prime divisor)
ignore ...

13432418150982799907
>>> n % p###line 4168:_sage_    : n % p
0
"""

def example_38():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> set_random_seed(Integer(0))###line 4203:_sage_    : set_random_seed(0)
>>> p = next_prime(randrange(Integer(2)**Integer(96)))###line 4204:_sage_    : p = next_prime(randrange(2^96))
>>> q = next_prime(randrange(Integer(2)**Integer(97)))###line 4205:_sage_    : q = next_prime(randrange(2^97))
>>> n = p * q###line 4206:_sage_    : n = p * q
>>> qsieve(n)###line 4207:_sage_    : qsieve(n)
([6340271405786663791648052309,
  46102313108592180286398757159], '')
"""

def example_39():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> legendre_symbol(Integer(2),Integer(3))###line 4397:_sage_    : legendre_symbol(2,3)
-1
>>> legendre_symbol(Integer(1),Integer(3))###line 4399:_sage_    : legendre_symbol(1,3)
1
>>> legendre_symbol(Integer(3),Integer(5))###line 4401:_sage_    : legendre_symbol(3,5)
-1
>>> legendre_symbol(Mod(Integer(3),Integer(5)), Integer(5))###line 4403:_sage_    : legendre_symbol(Mod(3,5), 5)
-1
"""

def example_40():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> legendre_symbol(Integer(69),Integer(389))###line 4573:_sage_    : legendre_symbol(69,389)
1
"""

def example_41():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def kr(a, p):###line 4632:_sage_    : def kr(a, p):
...    if Mod(a,p)**((p-Integer(1))//Integer(2)) == Integer(1):
...       return Integer(1)
...    else:
...       return -Integer(1)
>>> for a in range(Integer(1),Integer(5)):###line 4637:_sage_    : for a in range(1,5):
...    print a, kr(a,Integer(5))
1 1
2 -1
3 -1
4 1
"""

def example_42():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> p = Integer(726377359)###line 4683:_sage_    : p = 726377359
>>> Mod(Integer(3), p)**((p-Integer(1))//Integer(2))###line 4684:_sage_    : Mod(3, p)^((p-1)//2)
726377358
"""

def example_43():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def gauss(a, p):###line 4776:_sage_    : def gauss(a, p):
...    # make the list of numbers reduced modulo p
...    v = [(n*a)%p for n in range(Integer(1), (p-Integer(1))//Integer(2) + Integer(1))]
...    # normalize them to be in the range -p/2 to p/2
...    v = [(x if (x < p/Integer(2)) else x - p) for x in v]
...    # sort and print the resulting numbers
...    v.sort()
...    print v
...    # count the number that are negative
...    num_neg = len([x for x in v if x < Integer(0)])
...    return (-Integer(1))**num_neg
>>> gauss(Integer(2), Integer(13))###line 4787:_sage_    : gauss(2, 13)
[-5, -3, -1, 2, 4, 6]
-1
>>> legendre_symbol(Integer(2),Integer(13))###line 4790:_sage_    : legendre_symbol(2,13)
-1
>>> gauss(Integer(4), Integer(13))###line 4792:_sage_    : gauss(4, 13)
[-6, -5, -2, -1, 3, 4]
1
>>> legendre_symbol(Integer(4),Integer(13))###line 4795:_sage_    : legendre_symbol(4,13)
1
>>> gauss(Integer(2),Integer(31))###line 4797:_sage_    : gauss(2,31)
[-15, -13, -11, -9, -7, -5, -3, -1, 2, 4, 6, 8, 10, 12, 14]
1
>>> legendre_symbol(Integer(2),Integer(31))###line 4800:_sage_    : legendre_symbol(2,31)
1
"""

def example_44():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> K = CyclotomicField(Integer(5),names=('zeta',)); (zeta,) = K._first_ngens(Integer(1))###line 5026:_sage_    : K.<zeta> = CyclotomicField(5)
>>> zeta**Integer(5)###line 5027:_sage_    : zeta^5
1
>>> Integer(1)/zeta###line 5029:_sage_    : 1/zeta
-zeta^3 - zeta^2 - zeta - 1
"""

def example_45():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def gauss_sum(a,p):###line 5050:_sage_    : def gauss_sum(a,p):
...   K = CyclotomicField(p,names=('zeta',)); (zeta,) = K._first_ngens(Integer(1))
...   return sum(legendre_symbol(n,p) * zeta**(a*n)
...              for n in range(Integer(1),p))
>>> g2 = gauss_sum(Integer(2),Integer(5)); g2###line 5054:_sage_    : g2 = gauss_sum(2,5); g2
2*zeta^3 + 2*zeta^2 + 1
>>> g2.complex_embedding()###line 5056:_sage_    : g2.complex_embedding()
-2.2360679775 + 3.33066907388e-16*I
>>> g2**Integer(2)###line 5058:_sage_    : g2^2
5

Here, $g_2$ is initially output as a polynomial in $\zeta_5$,
so there is no loss of precision.  The \code{complex\_embedding}
command shows some embedding of $g_2$ into the complex numbers,
which is only correct to about the first 15 digits.  Note
that $g_2^2 = 5$, so $g_2 = -\sqrt{5}$.

We compute a graphical representation of the Gauss sum $g_2$
as follows (see Figure~\ref{fig:gauss_sum}):

zeta = CDF(exp(2*pi*I/5))
v = [legendre_symbol(n,5) * zeta^(2*n) for n in range(1,5)]
S = sum([point(tuple(z), pointsize=100) for z in v])
show(S + point(tuple(sum(v)), pointsize=100, rgbcolor='red'))

\end{sg}

\begin{figure}
\includegraphics[width=\textwidth]{graphics/gauss_sum}
\caption{The red dot is the Gauss sum $g_2$ for $p=5$\label{fig:gauss_sum}}
\end{figure}

Figure~\ref{fig:gauss_sum} illustrates the Gauss sum $g_2$ for $p=5$.
The Gauss sum is obtained by adding the points on the unit circle,
with signs as indicated, to obtain the real number $-\sqrt{5}$.  This
suggests the following proposition, whose proof will require some
work.

\begin{proposition}[Gauss Sum]\label{prop:gauss_sum1}\iprop{Gauss sum}
For any~$a$ not divisible by~$p$,
$$
\ds g_a^2 = (-1)^{(p-1)/2}p.
$$
\end{proposition}

\begin{sg}
We illustrate using \sage that the proposition is
correct for $p=7$ and $p=13$:

>>> [gauss_sum(a, Integer(7))**Integer(2) for a in range(Integer(1),Integer(7))]###line 5101:_sage_    : [gauss_sum(a, 7)^2 for a in range(1,7)]
[-7, -7, -7, -7, -7, -7]
>>> [gauss_sum(a, Integer(13))**Integer(2) for a in range(Integer(1),Integer(13))]###line 5103:_sage_    : [gauss_sum(a, 13)^2 for a in range(1,13)]
[13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
"""

def example_46():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> S = PolynomialRing(GF(Integer(13)),names=('x',)); (x,) = S._first_ngens(Integer(1))###line 5354:_sage_    : S.<x> = PolynomialRing(GF(13))
>>> R = S.quotient(x**Integer(2) - Integer(3),names=('alpha',)); (alpha,) = R._first_ngens(Integer(1))###line 5355:_sage_    : R.<alpha> = S.quotient(x^2 - 3)
>>> (Integer(2)+Integer(3)*alpha)*(Integer(1)+Integer(2)*alpha)###line 5356:_sage_    : (2+3*alpha)*(1+2*alpha)
7*alpha + 7
"""

def example_47():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def find_sqrt(a, p):###line 5414:_sage_    : def find_sqrt(a, p):
...    assert (p-Integer(1))%Integer(4) == Integer(0)
...    assert legendre_symbol(a,p) == Integer(1)
...    S = PolynomialRing(GF(p),names=('x',)); (x,) = S._first_ngens(Integer(1))
...    R = S.quotient(x**Integer(2) - a,names=('alpha',)); (alpha,) = R._first_ngens(Integer(1))
...    while True:
...        z = GF(p).random_element()
...        w = (Integer(1) + z*alpha)**((p-Integer(1))//Integer(2))
...        (u, v) = (w[Integer(0)], w[Integer(1)])
...        if v != Integer(0): break
...    if (-u/v)**Integer(2) == a: return -u/v
...    if ((Integer(1)-u)/v)**Integer(2) == a: return (Integer(1)-u)/v
...    if ((-Integer(1)-u)/v)**Integer(2) == a: return (-Integer(1)-u)/v
...
>>> b = find_sqrt(Integer(3),Integer(13))###line 5428:_sage_    : b = find_sqrt(3,13)
>>> print "ignore this";  b                        # random: either 9 or 3###line 5429:_sage_    : b                        # random: either 9 or 3
ignore ...

9
>>> b**Integer(2)###line 5431:_sage_    : b^2
3
>>> b = find_sqrt(Integer(3),Integer(13))###line 5433:_sage_    : b = find_sqrt(3,13)
>>> print "ignore this";  b                        # see, it's random###line 5434:_sage_    : b                        # see, it's random
ignore ...

4
>>> print "ignore this";  find_sqrt(Integer(5),Integer(389))         # random: either 303 or 86###line 5436:_sage_    : find_sqrt(5,389)         # random: either 303 or 86
ignore ...

303
>>> print "ignore this";  find_sqrt(Integer(5),Integer(389))         # see, it's random###line 5438:_sage_    : find_sqrt(5,389)         # see, it's random
ignore ...

86
"""

def example_48():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> continued_fraction(Integer(17)/Integer(23))###line 5702:_sage_    : continued_fraction(17/23)
[0, 1, 2, 1, 5]
>>> continued_fraction(e)###line 5704:_sage_    : continued_fraction(e)
[2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1,
 12, 1, 1, 11]
"""

def example_49():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> continued_fraction(e, bits=Integer(20))###line 5712:_sage_    : continued_fraction(e, bits=20)
[2, 1, 2, 1, 1, 4, 1, 1, 6]
>>> continued_fraction(e, bits=Integer(30))###line 5714:_sage_    : continued_fraction(e, bits=30)
[2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1]
"""

def example_50():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> a = continued_fraction(Integer(17)/Integer(23)); a###line 5720:_sage_    : a = continued_fraction(17/23); a
[0, 1, 2, 1, 5]
>>> a.value()###line 5722:_sage_    : a.value()
17/23
>>> b = continued_fraction(Integer(6)/Integer(23)); b###line 5724:_sage_    : b = continued_fraction(6/23); b
[0, 3, 1, 5]
>>> a + b###line 5726:_sage_    : a + b
[1]
"""

def example_51():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> c = continued_fraction(pi,bits=Integer(33)); c###line 5778:_sage_    : c = continued_fraction(pi,bits=33); c
[3, 7, 15, 1, 292, 2]
>>> c.convergents()###line 5780:_sage_    : c.convergents()
[3, 22/7, 333/106, 355/113, 103993/33102, 208341/66317]
"""

def example_52():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> c = continued_fraction(pi); c###line 5839:_sage_    : c = continued_fraction(pi); c
[3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
>>> for n in range(-Integer(1), len(c)):###line 5841:_sage_    : for n in range(-1, len(c)):
...    print c.pn(n)*c.qn(n-Integer(1)) - c.qn(n)*c.pn(n-Integer(1)),
1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1
>>> for n in range(len(c)):###line 5844:_sage_    : for n in range(len(c)):
...    print c.pn(n)*c.qn(n-Integer(2)) - c.qn(n)*c.pn(n-Integer(2)),
3 -7 15 -1 292 -1 1 -1 2 -1 3 -1 14 -3
"""

def example_53():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> c = continued_fraction([Integer(1),Integer(2),Integer(3),Integer(4),Integer(5)])###line 5866:_sage_    : c = continued_fraction([1,2,3,4,5])
>>> c.convergents()###line 5867:_sage_    : c.convergents()
[1, 3/2, 10/7, 43/30, 225/157]
>>> [c.pn(n) for n in range(len(c))]###line 5869:_sage_    : [c.pn(n) for n in range(len(c))]
[1, 3, 10, 43, 225]
>>> [c.qn(n) for n in range(len(c))]###line 5871:_sage_    : [c.qn(n) for n in range(len(c))]
[1, 2, 7, 30, 157]
"""

def example_54():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> c = continued_fraction([Integer(1),Integer(1),Integer(1),Integer(1),Integer(1),Integer(1),Integer(1),Integer(1)])###line 5905:_sage_    : c = continued_fraction([1,1,1,1,1,1,1,1])
>>> v = [(i, c.pn(i)/c.qn(i)) for i in range(len(c))]###line 5906:_sage_    : v = [(i, c.pn(i)/c.qn(i)) for i in range(len(c))]
>>> P = point(v, rgbcolor=(Integer(0),Integer(0),Integer(1)), pointsize=Integer(40))###line 5907:_sage_    : P = point(v, rgbcolor=(0,0,1), pointsize=40)
>>> L = line(v, rgbcolor=(RealNumber('0.5'),RealNumber('0.5'),RealNumber('0.5')))###line 5908:_sage_    : L = line(v, rgbcolor=(0.5,0.5,0.5))
>>> L2 = line([(Integer(0),c.value()),(len(c)-Integer(1),c.value())],       thickness=RealNumber('0.5'), rgbcolor=(RealNumber('0.7'),Integer(0),Integer(0)))###line 5910:_sage_    : L2 = line([(0,c.value()),(len(c)-1,c.value())],       thickness=0.5, rgbcolor=(0.7,0,0))
>>> (L+L2+P).show(xmin=Integer(0),ymin=Integer(1))###line 5911:_sage_    : (L+L2+P).show(xmin=0,ymin=1)
"""

def example_55():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def cf(bits):###line 6062:_sage_    : def cf(bits):
...   x = (Integer(1) + sqrt(RealField(bits)(Integer(5)))) / Integer(2)
...   return continued_fraction(x)
>>> cf(Integer(10))###line 6065:_sage_    : cf(10)
[1, 1, 1, 1, 1, 1, 1, 3]
>>> cf(Integer(30))###line 6067:_sage_    : cf(30)
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 2]
>>> cf(Integer(50))###line 6070:_sage_    : cf(50)
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
"""

def example_56():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def cf_sqrt_d(d, bits):###line 6529:_sage_    : def cf_sqrt_d(d, bits):
...   x = sqrt(RealField(bits)(d))
...   return continued_fraction(x)
>>> cf_sqrt_d(Integer(389),Integer(50))###line 6532:_sage_    : cf_sqrt_d(389,50)
[19, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38]
>>> cf_sqrt_d(Integer(389),Integer(100))###line 6534:_sage_    : cf_sqrt_d(389,100)
[19, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38,
 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1,
 2, 1, 1]
"""

def example_57():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def newton_root(f, iterates=Integer(2), x0=Integer(0), prec=Integer(53)):###line 6837:_sage_    : def newton_root(f, iterates=2, x0=0, prec=53):
...    x = RealField(prec)(x0)
...    R = PolynomialRing(ZZ,'x')
...    f = R(f)
...    g = f.derivative()
...    for i in range(iterates):
...        x = x - f(x)/g(x)
...    return x


\par\noindent{}Next we run the Newton iteration, and compute
the continued fraction of the result:

>>> a = newton_root(Integer(3847)*x**Integer(2) - Integer(14808904)*x + Integer(36527265)); a###line 6851:_sage_    : a = newton_root(3847*x^2 - 14808904*x + 36527265); a
2.46815700480740
>>> cf = continued_fraction(a); cf###line 6853:_sage_    : cf = continued_fraction(a); cf
[2, 2, 7, 2, 1, 5, 1, 1, 1, 1, 1, 1, 103, 8, 1, 2, 3, 1, 1]


\par\noindent{}We truncate the continued fraction and compute
its value.

>>> c = cf[:Integer(12)]; c###line 6861:_sage_    : c = cf[:12]; c
[2, 2, 7, 2, 1, 5, 1, 1, 1, 1, 1, 1]
>>> c.value()###line 6863:_sage_    : c.value()
9495/3847
"""

def example_58():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def sum_of_two_squares_naive(n):###line 6899:_sage_    : def sum_of_two_squares_naive(n):
...    for i in range(int(sqrt(n))):
...        if is_square(n - i**Integer(2)):
...            return i, (Integer(n-i**Integer(2))).sqrt()
...    return "%s is not a sum of two squares"%n

We next use our function in a couple of cases.

>>> sum_of_two_squares_naive(Integer(23))###line 6908:_sage_    : sum_of_two_squares_naive(23)
'23 is not a sum of two squares'
>>> sum_of_two_squares_naive(Integer(389))###line 6910:_sage_    : sum_of_two_squares_naive(389)
(10, 17)
>>> sum_of_two_squares_naive(Integer(2007))###line 6912:_sage_    : sum_of_two_squares_naive(2007)
'2007 is not a sum of two squares'
>>> sum_of_two_squares_naive(Integer(2008))###line 6914:_sage_    : sum_of_two_squares_naive(2008)
'2008 is not a sum of two squares'
>>> sum_of_two_squares_naive(Integer(2009))###line 6916:_sage_    : sum_of_two_squares_naive(2009)
(28, 35)
>>> Integer(28)**Integer(2) + Integer(35)**Integer(2)###line 6918:_sage_    : 28^2 + 35^2
2009
>>> sum_of_two_squares_naive(Integer(2)*Integer(3)**Integer(4)*Integer(5)*Integer(7)**Integer(2)*Integer(13))###line 6920:_sage_    : sum_of_two_squares_naive(2*3^4*5*7^2*13)
(189, 693)

\end{sg}

\begin{definition}[Primitive]
A representation $n=x^2 + y^2$ is \defn{primitive}\index{primitive!representation}
if~$x$ and~$y$ are coprime.
\end{definition}

\begin{lemma}\label{lem:primitive}
If~$n$ is divisible by a prime~$p\con 3\pmod{4}$, then~$n$
has no primitive representations.
\end{lemma}
\begin{proof}
Suppose $n$ has a primitive representation, $n=x^2 + y^2$, and
let~$p$ be any prime factor of~$n$.  Then
$$
   p \mid x^2 + y^2\quad \text{ and }\quad \gcd(x,y)=1,
$$
so $p\nmid x$ and $p\nmid y$.
Since $\zmod{p}$ is a field, we may divide by $y^2$ in the equation
$x^2 + y^2 \con 0\pmod{p}$ to see that
$
  (x/y)^2 \con -1\pmod{p}.
$
Thus the Legendre symbol $\kr{-1}{p}$ equals $+1$.
However, by Proposition~\ref{prop:euler},
$$
   \kr{-1}{p} = (-1)^{(p-1)/2}
$$
so $\kr{-1}{p}=1$ if and only if $(p-1)/2$ is even, which is
to say $p\con 1\pmod{4}$.
\end{proof}

\begin{proof}[Proof of Theorem~\ref{thm:sumsquare}
  $\left(\Longrightarrow\right)$]
  Suppose that $p\con 3\pmod{4}$ is a prime, that
  $p^r\mid n$ but $p^{r+1}\nmid n$ with~$r$ odd, and that
  $n=x^2+y^2$.  Letting $d=\gcd(x,y)$, we have
  $$
  x = dx', \quad y = dy', \quad\text{ and }\quad n = d^2 n'
  $$
  with $\gcd(x',y')=1$ and
  $$
  (x')^2 + (y')^2 = n'.
  $$
  Because~$r$ is odd, $p\mid n'$, so Lemma~\ref{lem:primitive}
  implies that $\gcd(x',y')>1$, which is a contradiction.
\end{proof}

To prepare for our proof of the implication
$\left(\Longleftarrow\right)$ of Theorem~\ref{thm:sumsquare}, we
reduce the problem to the case when~$n$ is prime.  Write $n=n_1^2 n_2$,
where $n_2$ has no prime factors $p\con 3\pmod{4}$.  It suffices to
show that~$n_2$ is a sum of two squares, since
\begin{equation}\label{eqn:ssnormprod}
  (x_1^2 + y_1^2)(x_2^2+y_2^2) = (x_1x_2-y_1y_2)^2 + (x_1y_2+x_2y_1)^2,
\end{equation}
so a product of two numbers that are sums of two squares is also
a sum of two squares.
Since~$2=1^2+1^2$ is a sum of two squares, it suffices to show
that any prime $p\con 1\pmod{4}$ is a sum of two squares.

\begin{lemma}\label{lem:approx}
If $x\in\R$ and $n\in\N$, then there is a fraction $\ds\frac{a}{b}$
in lowest terms such that $0<b\leq n$ and
$$\left| x - \frac{a}{b} \right| \leq \frac{1}{b(n+1)}.$$
\end{lemma}
\begin{proof}
Consider the continued fraction\index{continued fraction}
$[a_0,a_1,\ldots]$ of~$x$.
By Corollary~\ref{cor:cfconv}, for each~$m$
$$
 \left| x - \frac{p_m}{q_m}\right|
  < \frac{1}{q_m \cdot q_{m+1}}.
$$
Since $q_{m+1}\geq q_m + 1$ and $q_0=1$,
either there exists an~$m$ such that $q_m\leq n < q_{m+1}$, or the
continued fraction\index{continued fraction} expansion of~$x$ is finite and $n$ is larger
than the denominator of the rational number~$x$, in which case
we take $\frac{a}{b}=x$ and are done.  In the first
case,
$$
  \left| x - \frac{p_m}{q_m}\right|
   < \frac{1}{q_m \cdot q_{m+1}}
      \leq \frac{1}{q_m \cdot (n+1)},$$
so $\ds\frac{a}{b} = \frac{p_m}{q_m}$ satisfies the conclusion of
the lemma.
\end{proof}

\begin{proof}[Proof of Theorem~\ref{thm:sumsquare} $\left(\Longleftarrow\right)$]
As discussed above, it suffices to prove that any prime
$p\con 1\pmod{4}$ is a sum of two squares. Since $p\con 1\pmod{4}$,
$$
  (-1)^{(p-1)/2} = 1,
$$
Proposition~\ref{prop:euler} implies that
$-1$ is a square modulo~$p$; i.e., there exists~$r\in\Z$ such
that $r^2\con -1\pmod{p}$.
Lemma~\ref{lem:approx}, with $n=\lfloor \sqrt{p}\rfloor$
and $x=-\frac{r}{p}$,
implies that there are integers $a, b$ such that
$0<b<\sqrt{p}$ and
$$
 \left| -\frac{r}{p} - \frac{a}{b}\right| \leq\frac{1}{b(n+1)} < \frac{1}{b\sqrt{p}}.
$$
Letting $c = rb + pa$, we have that
$$
  |c| < \frac{pb}{b\sqrt{p}} = \frac{p}{\sqrt{p}} = \sqrt{p}
$$
so
$$
   0 < b^2 + c^2 < 2p.
$$
But $c \con rb\pmod{p}$, so
$$
  b^2 + c^2 \con b^2 + r^2 b^2 \con b^2(1+r^2) \con 0\pmod{p}.
$$
Thus $b^2 + c^2 = p$.
\end{proof}

\begin{remark}
Our proof of Theorem~\ref{thm:sumsquare} leads to an efficient
algorithm to compute a representation of any $p\con 1\pmod{4}$
as a sum of two squares.
\end{remark}

\begin{sg}
We next use \sage and Theorem~\ref{thm:sumsquare} to
give an efficient algorithm for writing a prime $p\con 1\pmod{4}$
as a sum of two squares.  First we implement the algorithm
that comes out of the proof of the theorem.

>>> def sum_of_two_squares(p):###line 7055:_sage_    : def sum_of_two_squares(p):
...    p = Integer(p)
...    assert p%Integer(4) == Integer(1), "p must be 1 modulo 4"
...    r = Mod(-Integer(1),p).sqrt().lift()
...    v = continued_fraction(-r/p)
...    n = floor(sqrt(p))
...    for x in v.convergents():
...        c = r*x.denominator() + p*x.numerator()
...        if -n <= c and c <= n:
...            return (abs(x.denominator()),abs(c))

Next we use the algorithm to write the first $10$-digit
prime $\con 1\pmod{4}$ as a sum of two squares:

>>> p = next_prime(next_prime(Integer(10)**Integer(10)))###line 7070:_sage_    : p = next_prime(next_prime(10^10))
>>> sum_of_two_squares(p)###line 7071:_sage_    : sum_of_two_squares(p)
(55913, 82908)

The above calculation was essentially instantanoues.
If instead we use the naive algorithm from before,
it takes several seconds to write $p$ as a sum of two squares.

>>> sum_of_two_squares_naive(p)###line 7079:_sage_    : sum_of_two_squares_naive(p)
(55913, 82908)
"""

def example_59():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> E = EllipticCurve([-Integer(5), Integer(4)])###line 7248:_sage_    : E = EllipticCurve([-5, 4])
>>> E###line 7249:_sage_    : E
Elliptic Curve defined by y^2  = x^3 - 5*x + 4
over Rational Field
>>> P = E.plot(thickness=Integer(4),rgbcolor=(RealNumber('0.1'),RealNumber('0.7'),RealNumber('0.1')))###line 7252:_sage_    : P = E.plot(thickness=4,rgbcolor=(0.1,0.7,0.1))
>>> P.show(figsize=[Integer(4),Integer(6)])###line 7253:_sage_    : P.show(figsize=[4,6])
"""

def example_60():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> E = EllipticCurve(GF(Integer(37)), [Integer(1),Integer(0)])###line 7270:_sage_    : E = EllipticCurve(GF(37), [1,0])
>>> E###line 7271:_sage_    : E
Elliptic Curve defined by y^2  = x^3 + x over
Finite Field of size 37
>>> E.plot(pointsize=Integer(45))###line 7274:_sage_    : E.plot(pointsize=45)
"""

def example_61():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> E = EllipticCurve([-Integer(5),Integer(4)])###line 7375:_sage_    : E = EllipticCurve([-5,4])
>>> P = E([Integer(1),Integer(0)]); Q = E([Integer(0),Integer(2)])###line 7376:_sage_    : P = E([1,0]); Q = E([0,2])
>>> P + Q###line 7377:_sage_    : P + Q
(3 : 4 : 1)
>>> P + P###line 7379:_sage_    : P + P
(0 : 1 : 0)
>>> P + Q + Q + Q + Q###line 7381:_sage_    : P + Q + Q + Q + Q
(350497/351649 : 16920528/208527857 : 1)
"""

def example_62():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> R = QQ['x1, y1, x2, y2, x3, y3, a, b']; (x1, y1, x2, y2, x3, y3, a, b,) = R._first_ngens(Integer(8))###line 7466:_sage_    : R.<x1,y1,x2,y2,x3,y3,a,b> = QQ[]

\vspace{1ex}

\noindent{}We define the relations the $x_i$ will satisfy, and a
quotient ring $Q$ in which those relations are satisfied.  (Quotients
of polynomial rings are a generalization of the construction $\Z/n\Z$
that may be viewed as the quotient of the ring $\Z$ of integers by the
relation that sets $n$ to equal $0$.)

>>> rels = [y1**Integer(2) - (x1**Integer(3) + a*x1 + b),###line 7477:_sage_    : rels = [y1^2 - (x1^3 + a*x1 + b),
...           y2**Integer(2) - (x2**Integer(3) + a*x2 + b),
...           y3**Integer(2) - (x3**Integer(3) + a*x3 + b)]
...
>>> Q = R.quotient(rels)###line 7481:_sage_    : Q = R.quotient(rels)

\vspace{1ex}

\noindent{}We define the group operation, which assumes the points are distinct.

>>> def op(P1,P2):###line 7488:_sage_    : def op(P1,P2):
...       x1,y1 = P1;  x2,y2 = P2
...       lam = (y1 - y2)/(x1 - x2); nu  = y1 - lam*x1
...       x3 = lam**Integer(2) - x1 - x2; y3 = -lam*x3 - nu
...       return (x3, y3)

\vspace{1ex}

\noindent{}We define three points, add them together via $P_1 + (P_2 + P_3)$
and $(P_1 + (P_2 + P_3))$, and observe that the results are the same modulo the relations.

>>> P1 = (x1,y1); P2 = (x2,y2); P3 = (x3,y3)###line 7500:_sage_    : P1 = (x1,y1); P2 = (x2,y2); P3 = (x3,y3)
>>> Z = op(P1, op(P2,P3)); W = op(op(P1,P2),P3)###line 7501:_sage_    : Z = op(P1, op(P2,P3)); W = op(op(P1,P2),P3)
>>> (Q(Z[Integer(0)].numerator()*W[Integer(0)].denominator() -###line 7502:_sage_    : (Q(Z[0].numerator()*W[0].denominator() -
...          Z[Integer(0)].denominator()*W[Integer(0)].numerator())) == Integer(0)
True
>>> (Q(Z[Integer(1)].numerator()*W[Integer(1)].denominator() -###line 7505:_sage_    : (Q(Z[1].numerator()*W[1].denominator() -
...          Z[Integer(1)].denominator()*W[Integer(1)].numerator())) == Integer(0)
True
"""

def example_63():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def lcm_upto(B):###line 7589:_sage_    : def lcm_upto(B):
...       return prod([p**int(math.log(B)/math.log(p))
...                    for p in prime_range(B+Integer(1))])
>>> lcm_upto(Integer(10)**Integer(2))###line 7592:_sage_    : lcm_upto(10^2)
69720375229712477164533808935312303556800
>>> LCM((ellipsis_range(Integer(1),Ellipsis,Integer(10)**Integer(2))))###line 7594:_sage_    : LCM([1..10^2])
69720375229712477164533808935312303556800
"""

def example_64():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def pollard(N, B=Integer(10)**Integer(5), stop=Integer(10)):###line 7724:_sage_    : def pollard(N, B=10^5, stop=10):
...       m = prod([p**int(math.log(B)/math.log(p))
...                 for p in prime_range(B+Integer(1))])
...       for a in (ellipsis_range(Integer(2),Ellipsis,stop)):
...           x = (Mod(a,N)**m - Integer(1)).lift()
...           if x == Integer(0): continue
...           g = gcd(x, N)
...           if g != Integer(1) or g != N: return g
...       return Integer(1)
>>> pollard(Integer(5917),Integer(5))###line 7733:_sage_    : pollard(5917,5)
61
>>> pollard(Integer(779167),Integer(5))###line 7735:_sage_    : pollard(779167,5)
1
>>> pollard(Integer(779167),Integer(15))###line 7737:_sage_    : pollard(779167,15)
2003
>>> pollard(Integer(4331),Integer(7))###line 7739:_sage_    : pollard(4331,7)
1
>>> pollard(Integer(4331),Integer(5))###line 7741:_sage_    : pollard(4331,5)
61
>>> pollard(Integer(187), Integer(15), Integer(2))###line 7743:_sage_    : pollard(187, 15, 2)
1
>>> pollard(Integer(187), Integer(15))###line 7745:_sage_    : pollard(187, 15)
11
"""

def example_65():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def ecm(N, B=Integer(10)**Integer(3), trials=Integer(10)):###line 7864:_sage_    : def ecm(N, B=10^3, trials=10):
...       m = prod([p**int(math.log(B)/math.log(p))
...                 for p in prime_range(B+Integer(1))])
...       R = Integers(N)
...       # Make Sage think that R is a field:
...       R.is_field = lambda : True
...       for _ in range(trials):
...           while True:
...               a = R.random_element()
...               if gcd(Integer(4)*a.lift()**Integer(3) + Integer(27), N) == Integer(1): break
...           try:
...               m * EllipticCurve([a, Integer(1)])([Integer(0),Integer(1)])
...           except ZeroDivisionError, msg:
...               # msg: "Inverse of <int> does not exist"
...               return gcd(Integer(str(msg).split()[Integer(2)]), N)
...       return Integer(1)
>>> set_random_seed(Integer(2))###line 7880:_sage_    : set_random_seed(2)
>>> ecm(Integer(5959), B=Integer(20))###line 7881:_sage_    : ecm(5959, B=20)
101
>>> ecm(next_prime(Integer(10)**Integer(20))*next_prime(Integer(10)**Integer(7)), B=Integer(10)**Integer(3))###line 7883:_sage_    : ecm(next_prime(10^20)*next_prime(10^7), B=10^3)
10000019
"""

def example_66():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> p = Integer(785963102379428822376694789446897396207498568951)###line 8089:_sage_    : p = 785963102379428822376694789446897396207498568951
>>> E = EllipticCurve(GF(p),     [Integer(317689081251325503476317476413827693272746955927),###line 8091:_sage_    : E = EllipticCurve(GF(p),     [317689081251325503476317476413827693272746955927,
...     Integer(79052896607878758718120572025718535432100651934)])
>>> E.cardinality()###line 8093:_sage_    : E.cardinality()
785963102379428822376693024881714957612686157429
>>> E.cardinality().is_prime()###line 8095:_sage_    : E.cardinality().is_prime()
True
>>> B = E([###line 8097:_sage_    : B = E([
...    Integer(771507216262649826170648268565579889907769254176),
...    Integer(390157510246556628525279459266514995562533196655)])
>>> n=Integer(670805031139910513517527207693060456300217054473)###line 8100:_sage_    : n=670805031139910513517527207693060456300217054473
>>> r=Integer(70674630913457179596452846564371866229568459543)###line 8101:_sage_    : r=70674630913457179596452846564371866229568459543
>>> P = E([Integer(14489646124220757767),###line 8102:_sage_    : P = E([14489646124220757767,
...    Integer(669337780373284096274895136618194604469696830074)])
>>> encrypt = (r*B, P + r*(n*B))###line 8104:_sage_    : encrypt = (r*B, P + r*(n*B))
>>> encrypt[Integer(1)] - n*encrypt[Integer(0)] == P   # decrypting works###line 8105:_sage_    : encrypt[1] - n*encrypt[0] == P   # decrypting works
True
"""

def example_67():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> T = lambda v: EllipticCurve(v###line 8245:_sage_    : T = lambda v: EllipticCurve(v
...          ).torsion_subgroup().invariants()
>>> T([-Integer(5),Integer(4)])###line 8247:_sage_    : T([-5,4])
[2]
>>> T([-Integer(43),Integer(166)])###line 8249:_sage_    : T([-43,166])
[7]
>>> T([-Integer(4),Integer(0)])###line 8251:_sage_    : T([-4,0])
[2, 2]
>>> T([-Integer(1386747), Integer(368636886)])###line 8253:_sage_    : T([-1386747, 368636886])
[8, 2]
"""

def example_68():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> r = lambda v: EllipticCurve(v).rank()###line 8270:_sage_    : r = lambda v: EllipticCurve(v).rank()
>>> r([-Integer(5),Integer(4)])###line 8271:_sage_    : r([-5,4])
1
>>> r([Integer(0),Integer(1)])###line 8273:_sage_    : r([0,1])
0
>>> r([-Integer(3024), Integer(46224)])###line 8275:_sage_    : r([-3024, 46224])
2
>>> r([-Integer(112), Integer(400)])###line 8277:_sage_    : r([-112, 400])
3
>>> r([-Integer(102627), Integer(12560670)])###line 8279:_sage_    : r([-102627, 12560670])
4
"""

def example_69():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


>>> def cong(n):###line 8472:_sage_    : def cong(n):
...       G = EllipticCurve([-n**Integer(2),Integer(0)]).gens()
...       if len(G) == Integer(0): return False
...       x,y,_ = G[Integer(0)]
...       return ((n**Integer(2)-x**Integer(2))/y,-Integer(2)*x*n/y,(n**Integer(2)+x**Integer(2))/y)
>>> cong(Integer(6))###line 8477:_sage_    : cong(6)
(3, 4, 5)
>>> cong(Integer(5))###line 8479:_sage_    : cong(5)
(3/2, 20/3, 41/6)
>>> cong(Integer(1))###line 8481:_sage_    : cong(1)
False
>>> cong(Integer(13))###line 8483:_sage_    : cong(13)
(323/30, 780/323, 106921/9690)
>>> (Integer(323)/Integer(30) * Integer(780)/Integer(323))/Integer(2)###line 8485:_sage_    : (323/30 * 780/323)/2
13
>>> (Integer(323)/Integer(30))**Integer(2) + (Integer(780)/Integer(323))**Integer(2) == (Integer(106921)/Integer(9690))**Integer(2)###line 8487:_sage_    : (323/30)^2 + (780/323)^2 == (106921/9690)^2
True
"""


if __name__ ==  '__main__':
    import doctest, sys
    s = doctest.testmod(sys.modules[__name__],
                   optionflags=doctest.NORMALIZE_WHITESPACE
                              |doctest.ELLIPSIS,
                   verbose=False,
                   globs=globals())
    quit_sage(verbose=False)
    sys.exit(s[0])
