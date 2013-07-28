class SavingDict(collections.MutableMapping):
    def __init__(self,
            filename,
            values      = None,
            default     = _raise,
            paired      = None):
        r"""
        dictionary-like class that saves itself on the filesystem

        INPUT:

        - ``filename`` -- file to store SavingDict

        - ``values`` -- dictionary-like object that sets initial values

          by default initial values will be loaded from filename, if it
          exists, otherwise the initial values will be empty

        - ``default`` -- callabable object requiring no arguments that
          provides a default value for keys not in SavingDict

        - ``paired`` -- another SavingDict that will be paired as per
          :meth:`set_paired`

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: sd = SavingDict(tmp)
            sage: sd['cat'] = 'meow'; sd['cat']
            'meow'
            sage: sd['cow'] = 'moo'; sd
            {'cow': 'moo', 'cat': 'meow'}
            sage: del sd; sd
            Traceback (most recent call last):
            ...
            NameError: name 'sd' is not defined
            sage: sd = SavingDict(tmp); sd
            {'cow': 'moo', 'cat': 'meow'}
            sage: os.unlink(tmp)
        """
        if not callable(default):
            raise ValueError("default must be callable")

        self._filename = filename
        if values is None:
            self._dict = load_dict_from_file(filename)
        else:
            self._dict = dict(values) # explicitly make copy
        self._default = default
        self._pairing = None
        if paired is not None:
            self.set_paired(paired)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: sd = SavingDict('', {0:1, 1:2})
            sage: repr(sd)
            '{0: 1, 1: 2}'
        """
        return repr(self._dict)

    def unset_pairing(self):
        r"""
        unset any pairing that was constructed

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp1, tmp2 = tempfile.mkstemp()[1], tempfile.mkstemp()[1]
            sage: sd1= SavingDict(tmp1); sd2 = SavingDict(tmp2, paired=sd1)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: sd1.unset_pairing()
            sage: del sd1[0]; sd1[0]
            Traceback (most recent call last):
            ...
            KeyError: 0
            sage: sd2[1]
            0
            sage: os.unlink(tmp1); os.unlink(tmp2)
        """
        if self._pairing:
            with self._pairing as paired:
                paired.unset_pairing()
        self._pairing = None

    def set_paired(self, other):
        r"""
        set another SavingDict to be updated with the reverse of this one
        and vice versa

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp1, tmp2 = tempfile.mkstemp()[1], tempfile.mkstemp()[1]
            sage: sd1, sd2 = SavingDict(tmp1), SavingDict(tmp2)
            sage: sd1.set_paired(sd2)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: del sd1[0]; sd1[0]
            Traceback (most recent call last):
            ...
            KeyError: 0
            sage: sd2[1]
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: sd2[2] = 3; sd1[3]; sd2[2]
            2
            3
            sage: sd2[2] = 4; sd1[3]
            Traceback (most recent call last):
            ...
            KeyError: 3
            sage: sd1[4]; sd2[2]
            2
            4
            sage: os.unlink(tmp1); os.unlink(tmp2)
        """
        if not isinstance(other, SavingDict):
            raise ValueError("other is not a SavingDict")

        self.unset_pairing()
        other.unset_pairing()

        class Pairing(object):
            def __init__(this, obj):
                this._obj = obj
            def __enter__(this):
                this._entered[0] += 1
                return other if this._obj is self else self
            def __exit__(this, type, value, tb):
                this._entered[0] -= 1
            def __nonzero__(this):
                return not this._entered[0]
        Pairing._entered = [0] # ints are immutable

        self._pairing, other._pairing = Pairing(self), Pairing(other)

    def _write(self):
        r"""
        writes self to disk

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: sd = SavingDict(tmp, {0:1, 1:2})
            sage: SavingDict(tmp)
            {}
            sage: sd._write()
            sage: SavingDict(tmp)
            {0: 1, 1: 2}
            sage: os.unlink(tmp)
        """
        tmpfile = self._filename + '%016x'%(random.randrange(256**8))
        s = cPickle.dumps(self._dict, protocol=2)
        with open(tmpfile, 'wb') as F:
            F.write(s)
        # This move is atomic (the files are on the same filesystem)
        os.rename(tmpfile, self._filename)

    def __setitem__(self, key, value):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: sd = SavingDict(tmp)
            sage: sd['cow'] = 'moo'
            sage: sd
            {'cow': 'moo'}
            sage: os.unlink(tmp)
        """
        if self._pairing:
            with self._pairing as paired:
                if key in self:
                    del paired[self[key]]
                paired[value] = key
        self._dict[key] = value
        self._write()

    def __delitem__(self, key):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: sd = SavingDict(tmp, {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: del sd['cow']
            sage: del sd['cow']
            Traceback (most recent call last):
            ...
            KeyError: 'cow'
            sage: sd
            {}
            sage: os.unlink(tmp)
        """
        if self._pairing and key in self:
            with self._pairing as paired:
                del paired[self[key]]
        del self._dict[key]
        self._write()

    def __getitem__(self, key):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: sd = SavingDict('', {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: sd['cow']
            'moo'
            sage: sd['moo']
            Traceback (most recent call last):
            ...
            KeyError: 'moo'
        """
        try:
            return self._dict[key]
        except KeyError:
            return self._default()

    def __contains__(self, key):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: sd = SavingDict('', {'cow': 'moo'}); sd
            {'cow': 'moo'}
            sage: 'cow' in sd
            True
            sage: 'moo' in sd
            False
        """
        return key in self._dict

    def __len__(self):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp = tempfile.mkstemp()[1]
            sage: sd = SavingDict(tmp); len(sd)
            0
            sage: sd['cow'] = 'moo'
            sage: len(sd)
            1
            sage: os.unlink(tmp)
        """
        return len(self._dict)

    def __iter__(self):
        r"""
        TESTS::

            sage: from sage.dev.git_interface import SavingDict
            sage: sd = SavingDict('', {'cow':'moo', 0:1}); sd
            {0: 1, 'cow': 'moo'}
            sage: for key in sd:
            ...       print key, sd[key]
            0 1
            cow moo
        """
        return iter(self._dict)

