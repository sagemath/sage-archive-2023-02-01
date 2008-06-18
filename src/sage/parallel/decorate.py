
class parallel:
    def __init__(self,
                 peval = 'reference',
                 threads = 2,
                 db = None,
                 compress = True):
        if isinstance(peval, str):
            if peval == 'reference':
                from reference import parallel_eval
                self.peval = parallel_eval
            else:
                raise ValueError, "Unknown peval system '%s'"%peval
        else:
            self.peval = peval
        self.threads = threads
        self.db = db
        self.compress = compress

    def __call__(self, f):
        return do_parallel(f = f, db = self.db, compress = self.compress,
                           threads = self.threads, peval = self.peval)

def do_parallel(f, db, compress, threads, peval):
    def g(*args, **kwds):
        if len(args) > 0 and isinstance(args[0], list):
            return peval(f, args[0], db=db, compress=compress, threads=threads)
        else:
            return f(*args, **kwds)
    return g
