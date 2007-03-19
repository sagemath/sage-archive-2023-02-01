"""
TESTS:

sage: R = Zp(5, prec=5, type='fixed-mod')
sage: a = random_matrix(R,5)
sage: a.determinant()                # random output
5 + 3*5^2 + 5^3 + 4*5^4 + O(5^5)
sage: K = Qp(3, 10,'capped-rel'); K.krull_dimension()
0


"""
