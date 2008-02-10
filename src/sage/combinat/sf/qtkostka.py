import sage.combinat.partition


def qt_kostka(n):
    K = {}

    sp = sage.combinat.partition.Partition_n(n).list()
    spc = map(lambda x: x.conjugate(), sp)

    for i in range(len(sp)):
        mu = spc[i]
        nu = sp[i]
        j = spc.index(nu)
        if j < i:
            for k in range(len(sp)):
                K[spc[k],mu] = K[sp[k],nu].subs(q=t,t=q)
        else:
            pass


def e_to_s(part):
    nu = sage.combinat.partition.Partition_class(part).conjugate()

