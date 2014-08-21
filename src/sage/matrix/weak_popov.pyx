from collections import defaultdict

# desc: gives last position of highest degree of polynomial in vector (so called "Leading position")
cdef leading_position(v):
	p = -1 # pos of max
	m = -1 # max
	for c in range(v.degree()):
		if(v[c].degree()>=m):
			m = v[c].degree()
			p = c
	return p


# desc: makes a simple transformation of M[rp2] as (will be) defined in documentation of mathpart (i=rp1,j=rp2)
# anot: rp2 is the one to be changed!
# cond: deg rp2 >= deg rp1
cdef simple_transformation(M,rp2,rp1,LP):
	delta = M[rp2][LP].degree()-M[rp1][LP].degree()
	alpha = (M[rp2][LP].coefficients()[-1]) / (M[rp1][LP].coefficients()[-1])
	for i in range(M.ncols()):
		M[rp2,i] -= alpha*M[rp1,i].shift(delta)
	#M[rp2] -= alpha*M[rp1].shift(delta)
	return

# desc: transforms M into weak-popov form
cpdef mulders_storjohann(M):
	lps = defaultdict(list)
	for c in range(M.nrows()):
		lp = leading_position(M[c])
		if not M[c,lp]==-1
			lps[lp].append(c)
	
	while len(lps)<M.nrows():
		for pos in lps:
			if len(lps[pos])>1:
				if (M[lps[pos][0]][pos].degree() >= M[lps[pos][1]][pos].degree()):
					arownr = lps[pos][0]
					brownr = lps[pos][1]
				else:
					arownr = lps[pos][1]
					brownr = lps[pos][0]
				simple_transformation(M,arownr,brownr,pos)
				lps[pos].remove(arownr)
				lps[leading_position(M[arownr])].append(arownr) 
				break
	return M
    
