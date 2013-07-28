from sage.env import DOT_SAGE

class Config(collections.MutableMapping):
    r"""
    Wrapper around the ``devrc`` file storing the configuration for
    :class:`SageDev`.

    INPUT:

    - ``devrc`` -- a string (default: the absolute path of the ``devrc`` file in ``DOT_SAGE``)

    EXAMPLES::

        sage: dev._config
        Config('''
        [trac]
        username = doctest
        password_timeout = .5
        ''')
    """
    def __init__(self, devrc = os.path.join(DOT_SAGE, 'devrc')):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.sagedev import Config
            sage: type(Config())
        """
        self._config = configparser.ConfigParser()
        self._devrc = devrc
        self._read_config()
        self._doctest_config = False

    def __repr__(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: repr(dev._config)
            "Config('''\n[trac]\nusername = doctest\npassword_timeout = .5\n''')"
        """
        return "Config('''\n"+"\n".join([ "[%s]\n"%s+"\n".join(["%s = %s"%(o,self[s][o]) for o in self[s] ]) for s in self ])+"\n''')"

    def _read_config(self):
        r"""
        Read the configuration from disk.

        EXAMPLES::

            sage: from sage.dev.sagedev import Config, doctest_config
            sage: c = doctest_config()
            sage: c._write_config()
            sage: c = Config(c._devrc)
            sage: c._read_config()
            sage: c
            Config('''
            [trac]
            username = doctest
            password_timeout = .5
            ''')
        """
        if os.path.exists(self._devrc):
            self._config.read(self._devrc)

    def _write_config(self):
        r"""
        Write the configuration to disk.

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: os.unlink(c._devrc)
            sage: os.path.exists(c._devrc)
            False
            sage: c._write_config()
            sage: os.path.exists(c._devrc)
            True
        """
        with open(self._devrc, 'w') as F:
            self._config.write(F)
        # set the configuration file to read only by this user,
        # because it may contain the trac password
        os.chmod(self._devrc, 0600)

    def __getitem__(self, section):
        r"""
        Return the configurations in ``section``.

        EXAMPLES::

            sage: dev._config['trac']
            IndexableForSection('''
            username = doctest
            password_timeout = .5
            ''')
            sage: dev._config['tig']
            Traceback (most recent call last):
            ...
            KeyError: 'tig'
        """
        if not section in self:
            raise KeyError(section)

        class IndexableForSection(collections.MutableMapping):
            def __init__(this, section):
                this._section = section
            def __repr__(this):
                return "IndexableForSection('''\n"+"\n".join(["%s = %s"%(o,this[o]) for o in this])+"\n''')"
            def __getitem__(this, option):
                try:
                    return self._config.get(this._section, option)
                except configparser.NoOptionError:
                    raise KeyError(option)
            def __iter__(this):
                return iter(self._config.options(this._section))
            def __setitem__(this, option, value):
                self._config.set(this._section, option, value)
                self._write_config()
            def getboolean(this, option):
                return self._config.getboolean(this._section, option)
            def __delitem__(this, option):
                self._config.remove_option(this._section, option)
                self._write_config()
            def __len__(this):
                return len(self._config.options(this._section))
            def __contains__(this, option):
                return option in self._config.options(this._section)

        return IndexableForSection(section)

    def __contains__(self, section):
        r"""
        returns true if section is in the configuration

        EXAMPLES::

            sage: 'trac' in dev._config
            True
            sage: 'nottrac' in dev._config
            False
        """
        return section in self._config.sections()

    def __iter__(self):
        r"""
        Return an iterator over the section names.

        EXAMPLES::

            sage: list(dev._config)
            ['trac']
        """
        return iter(self._config.sections())

    def __setitem__(self, section, dictionary):
        r"""
        Set ``section`` to ``dictionary``.

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: c['foo'] = {'foo':'foo'}
            sage: c['foo']['foo']
            'foo'
        """
        if self._config.has_section(section):
            self.remove_section(section)
        self._config.add_section(section)
        for option, value in dictionary.iteritems():
            self._config.set(section, option, value)
        self._write_config()

    def __len__(self):
        r"""
        get the number of sections in the configuration

        EXAMPLES::

            sage: len(dev._config)
            1
        """
        return len(self._config.sections())

    def __delitem__(self, section):
        r"""
        remove ``section`` from the configuration

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: del c['git']
            sage: list(c)
            ['trac']
        """
        self._config.remove_section(section)
        self._write_config()
