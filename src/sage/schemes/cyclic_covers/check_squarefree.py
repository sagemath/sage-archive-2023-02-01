def check_squarefree(f):
    should_be_coprime = [f, f.derivative()]
    try:
        squarefree = should_be_coprime[0].gcd(should_be_coprime[1]).degree()==0
    except (AttributeError, NotImplementedError, TypeError):
        try:
            squarefree = should_be_coprime[0].resultant(should_be_coprime[1])!=0
        except (AttributeError, NotImplementedError, TypeError):
            raise NotImplementedError("Cannot determine whether " \
                      "polynomial %s is square free." % (f,));
    return squarefree
