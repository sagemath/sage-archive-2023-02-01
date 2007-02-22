"""nodoctest
Generic database that uses ZODB.
"""



"""
Important long email note.  This guy found that zlib is much better
than bzip2 to compression of small files with ZODB.

Hi William,

than you might be interested in what I found in the meantime, when doing some
simple comparisons for a particular application (ran it as "time
my-application.py"):

    tests on local disk

    no compression
    real    28m6.191s
    user    24m44.077s
    sys     0m54.764s
    resulting Data.fs: 1.9G

    ----------

    bzip2, thresh 1024
    real    53m36.815s
    user    49m1.959s
    sys     1m6.094s
    resulting Data.fs: 871M

    ----------

    bzip2, thresh 2048
    real    55m52.140s
    user    49m37.324s
    sys     1m6.361s
    resulting Data.fs: 871M

    ----------

    zlib, thresh 1024
    real    34m44.240s
    user    30m7.113s
    sys     0m42.538s
    resulting Data.fs: 852M

    ----------

    zlib, thresh 2048
    real    32m38.335s
    user    30m21.959s
    sys     0m42.355s
    resulting Data.fs: 852M


I found that to be very interesting: in my case, zlib compresses a little bit
better than bzip2 (interesting enough, but I think bzip2 has its strength on
larger chunks anyway), but more important for me, seems to work nearly as
fast as the plain FileStorage. The test was not run under laboratory
conditions, but should still be good enough to get a trend.

Cheers,

Sascha

--
Gallileus - the power of knowledge

Gallileus GmbH                   http://www.gallileus.info/

"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os, shutil, cPickle

import transaction
import BTrees.OOBTree
from ZODB import FileStorage, DB

import sage.databases.compressed_storage
import sage.misc.misc

# The following effectively turns of the ZODB logger, which is OK for us.
# Without this, one gets this annoying error message a lot:
#       No handlers could be found for logger "ZODB.FileStorage"
import logging
logging.getLogger("ZODB.FileStorage").setLevel(10000000)
logging.getLogger("ZODB.lock_file").setLevel(10000000)
logging.getLogger("ZODB.Connection").setLevel(10000000)

#import Globals, ZODB.FileStorage, ZODB.compressed_storage

Storage = sage.databases.compressed_storage.CompressedStorage(FileStorage.FileStorage, 1024)

DB_HOME = "%s/data/"%sage.misc.misc.SAGE_ROOT

class _uniq(object):
    _db = {} # Class variable, no globals!

    def __new__(cls, name="", read_only=True, unique_key=None):
        key = (cls, unique_key)
        if _uniq._db.has_key(key):
            return _uniq._db[key]
        X = object.__new__(cls)
        _uniq._db[key] = X
        return X

class Database(_uniq):
    def __init__(self, name, read_only=True, thresh=1024):
        if not hasattr(self, 'name'):
            self.read_only = read_only
            self.name = name
            self._thresh = thresh
            self._load_()

    def _load_(self):
        name = self.name
        read_only = self.read_only
        thresh = self._thresh
        if not os.path.exists("%s/%s"%(DB_HOME,name)):
            try:
                os.makedirs("%s/%s"%(DB_HOME,name))
            except OSError:    # for online calculator...
                pass
        self._dbname = "%s/%s/%s"%(DB_HOME, name, name)
        if self.read_only and not os.path.exists(self._dbname):
            raise RuntimeError, "The database %s is not installed."%self._dbname
        fs = FileStorage.FileStorage(self._dbname, read_only=self.read_only)
        self._storage = sage.databases.compressed_storage.CompressedStorage(fs, thresh=self._thresh)
        self._db = DB(self._storage)
        self.conn = self._db.open()
        self._root = self.conn.root()
        if not self._root.has_key("btree"):
            self._root["btree"] = BTrees.OOBTree.OOBTree()
        self.root = self._root["btree"]

    def begin(self):
        r"""Start a new database transaction"""
        transaction.get().begin()

    def abort(self):
        r"""Abort the current database transaction, without committing"""
        transaction.get().abort()

    def commit(self):
        """
        Commit the new version of this object to the database file.

        Note that if a data item corresponding to a key is changed,
        you still have to tell the database that that data item
        was changed by calling the changed method with that key.
        """
        if self.read_only:
            raise RuntimeError, "Cannot commit read only database."
        self._root._p_changed = 1
        transaction.get().commit()
        #get_transaction().commit()

    def changed(self, key):
        """
        Informs the database that some items corresponding to
        the given key may have changed.  This does not commit
        the changes to disk (use the commit function after
        calling changed to do that).
        """
        self.root._p_changed = 1
        X = self.root[key]
        self.root[key] = X

    def pack(self):
        """
        This function is not implemented -- I couldn't get pack
        working with compressed storage.  You can use the rebuild
        function instead, though it's slower than the usual ZODB pack,
        since it completely rebuilds the database from scratch.
        """
        raise NotImplementedError
        self._db.pack()
        self.commit()

    def rebuild(self, thresh=None):
        """
        Completely rebuild the database from scratch, by going
        through and writing everything out to a temporary database,
        then moving the temporary database files over self's
        files.  This can take a long time.

        The main reason for this function is that unfortunately I
        can't get pack to work on compressed ZODB databases.

        A copy of the old database file is created before rebuild.

        If you specify a thresh then that threshhold is used for
        recompressing all the objects.  Note that the threshhold is
        not saved as part of the database, so new objects will be
        compressed using whatever threshhold you use when creating
        the database object.
        """
        if self.read_only:
            raise RuntimeError, "Cannot pack read only database."
        if thresh == None:
            thresh = self._thresh
        else:
            self._thresh = thresh
        rebuild_name = self._dbname + "_rebuild"
        shutil.copy2(self._dbname, self._dbname + ".old")
        if os.path.exists(rebuild_name):
            os.unlink(rebuild_name)
        fs = FileStorage.FileStorage(rebuild_name, read_only=False)
        storage = sage.databases.compressed_storage.CompressedStorage(fs, thresh)
        db = DB(storage)
        conn = db.open()
        _root = conn.root()
        root = BTrees.OOBTree.OOBTree()
        _root["btree"] = root
        for k, x in self.root.iteritems():
            root[k] = x
        _root._p_changed = 1
        #get_transaction().commit()
        transaction.get().commit()
        shutil.move(rebuild_name, self._dbname)
        os.unlink(rebuild_name + ".tmp")
        os.unlink(rebuild_name + ".index")
        os.unlink(rebuild_name + ".lock")
        self.read_only = True


    def __repr__(self):
        return "Database %s"%self.name

    def __setitem__(self, x, y):
        try:
            self.root[x] = y
        except AttributeError:
            self._init()
            self.root[x] = y

    def __getitem__(self, x):
        try:
            if not isinstance(x, slice):
                return self.root[x]
            return [self[k] for k in range(x.start, x.stop, x.step)]
        except AttributeError:
            self._init()
            return self.root[x]

    def __delitem__(self, x):
        del self.root[x]

    def has_key(self, x):
        return bool(self.root.has_key(x))

    def keys(self):
        return self.root.keys()

    def as_dict(self, keys=None):
        """
        Return a dict representation of the database.

        Since the database could be large, if the optional keys
        parameter is given then only the elements of the database
        with key in keys are listed.
        """
        X = {}
        if keys == None:
            keys = self.root.keys()
        for k in keys:
            if self.has_key(k):
                X[k] = self.root[k]
        return X

    def dump_as_dict(self, filename, keys):
        X = self.as_dict(keys)
        print "Dumping %s..."%filename
        s = cPickle.dumps(X,2)
        dir = "%s/pickles/"%DB_HOME
        if not os.path.exists(dir):
            os.makedirs(dir)
        open("%s/%s"%(dir,filename), "w").write(s)

    def dump_as_dict_intervals(self, basename, Nstart, Nstop, length):
        N = Nstart
        while N <= Nstop:
            N2 = min(Nstop, N+length)
            Z = xrange(N, N2+1)
            self.dump_as_dict("%s_%s-%s"%(basename,N,N2), Z)
            N += length

    def restore_from_dict(self, filename):
        """
        Restore from the filename which must store a pickled dict.

        After loading the database is committed.
        """
        if self.read_only:
            raise RuntimeError, "%s is read only."%self
        dir = "%s/pickles/"%DB_HOME
        s = open("%s/%s"%(dir,filename)).read()
        print "Restoring %s..."%filename
        X = cPickle.loads(s,2)
        for k, x in X.iteritems():
            self.root[k] = x
        self.commit()

    def restore_from_dict_all(self, basename):
        """
        Restore all files that start with the given basename.

        Each file is loaded then commited to disk before the next
        file is loaded.
        """
        X = os.listdir("%s/pickles/"%DB_HOME)
        n = len(basename)
        for F in X:
            if F[:n] == basename:
                self.restore_from_dict(F)

    def delete_all(self):
        """
        Delete every entry in the database.
        """
        del self._root["btree"]
        self._root["btree"] = BTree.BTree()
        self.root = self._root["btree"]

    def clone(self, new_name):
        """
        Copy the database to a new database with the given new_name.
        There must not be a database with the new_name already, or a
        RuntimeError exception is raised.
        """
        if os.path.exists("%s/%s"%(DB_HOME,new_name)):
            raise RuntimeError, "Cannot clone to %s since that database already exists."%name
        os.path.makedirs("%s/%s"%(DB_HOME,new_name))
        shutil.copy2("%s/%s/%s"%(DB_HOME,name,name), "%s/%s"%(DB_HOME,new_name))
