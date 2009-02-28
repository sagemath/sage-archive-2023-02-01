
import bz2, os

import sage.misc.misc

import sage.databases.db   # very important that this be fully qualified


DB_HOME = "%s/data/kohel/"%sage.misc.misc.SAGE_ROOT


def dbz_to_intlist(name):
    file = DB_HOME + name
    if not os.path.exists(file):
        raise RuntimeError, "Kohel database file %s not available"%name
    data = bz2.decompress(open(file).read())
    # this line assumes a newline at end of file
    data = "[[" + data.replace("\n","],[").replace(" ",",")[:-2] + "]"
    return eval(data)


def _padzeros(s):
    return "0"*(6-len(str(s))) + str(s)

class BrandtModule:
    def __init__(self, level):
        start = 400*(level//401) + 1
        stop = start + 399

        base = "ModBrdt/%s-%s/%s/"%(_padzeros(start), _padzeros(stop),
                                    _padzeros(level))
        padded_level = _padzeros(level)
        brdt = "%s/%s.%s.001.brdt.dbz"%(base,padded_level,padded_level)
        self.level = level
        try:
            self.raw = dbz_to_intlist(brdt)
        except RuntimeError:
            raise RuntimeError, "No data in the database about BrandtModule of level %s"%level


        self.auts = []
        tmp = self.raw[0]
        for i in range(len(tmp)/2):
            self.auts += [tmp[2*i+1]] * tmp[2*i]

        self.dimension = len(self.auts)


    def __repr__(self):
        return "Brandt Module of Level %s of Dimension %s"%(\
            self.level, self.dimension)


class KohelBrandtModuleDatabase(sage.databases.db.Database):
    def __init__(self, read_only=True):
        """
        Initialize the database.

        INPUT:
            read_only -- bool (default: True), if True, then the
                         database is read_only and changes cannot be
                         committed to disk.
        """
        sage.databases.db.Database.__init__(self,
             name="kohel_brandt_modules", read_only=read_only)

    def __repr__(self):
        return "Kohel database"

    def load(self, level):
        if self.read_only:
            raise RuntimeError, "Cannot load new module into read only database."
        B = BrandtModule(level)
        self[level] = B
        self.commit()


