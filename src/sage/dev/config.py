r"""
Configuration wrapper

A wrapper around the ``devrc`` file in ``DOT_SAGE`` directory which stores the
configuration for :class:`SageDev`.

AUTHORS:

- TODO: add authors from github's history and trac's history

"""
#*****************************************************************************
#       Copyright (C) 2013 TODO
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.env import DOT_SAGE

class Config(collections.MutableMapping):
    r"""
    Wrapper around the ``devrc`` file storing the configuration for
    :class:`SageDev`.

    INPUT:

    - ``devrc`` -- a string (default: the absolute path of the ``devrc`` file
    in ``DOT_SAGE``)

    EXAMPLES::

        sage: from sage.dev.sagedev import Config
        sage: Config() # random output, the output depends on the user's config
        sage: dev._config
        Config('''
        [trac]
        username = doctest
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

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: c["trac"]["username"] = "foo"

            sage: from sage.dev.sagedev.config import Config
            sage: c2 = Config(c._devrc)
            sage: c2
            sage: c._write_config()
            sage: c2._read_config()
            sage: c2
            Config('''
            [trac]
            username = foo
            ''')

        """
        if os.path.exists(self._devrc):
            self._config.read(self._devrc)

    def _write_config(self):
        r"""
        Write the configuration to disk.

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
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

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: c['trac']
            IndexableForSection('''
            username = doctest
            ''')
            sage: c['tig']
            Traceback (most recent call last):
            ...
            KeyError: 'tig'

        """
        if not section in self:
            raise KeyError(section)

        class IndexableForSection(collections.MutableMapping):
            r"""
            A section of a :class:`Config`.

            TESTS::

                sage: from sage.dev.sagedev.test.config import Config
                sage: c = Config()
                sage: c["trac"]

            """
            def __init__(this, section):
                r"""
                Initialization.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: type(c['trac'])

                """
                this._section = section

            def __repr__(this):
                r"""
                Return a printable representation of this section.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: c["trac"]

                """
                return "IndexableForSection('''\n"+"\n".join(["%s = %s"%(o,this[o]) for o in this])+"\n''')"

            def __getitem__(this, option):
                r"""
                Return the value for ``option`` or raise a ``KeyError`` if
                this option is not set.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: c["trac"]["username"]
                    doctest
                    sage: c["trac"]["nousername"]
                    Traceback (most recent call last):
                    ...
                    KeyError: nousername

                """
                try:
                    return self._config.get(this._section, option)
                except configparser.NoOptionError:
                    raise KeyError(option)

            def __iter__(this):
                r"""
                Return an iterator over the options in this section.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: list(c["trac"])

                """
                return iter(self._config.options(this._section))

            def __setitem__(this, option, value):
                r"""
                Set ``option`` to ``value``.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: c["trac"]["username"] = "foo"
                    sage: c["trac"]["username"]
                    foo

                """
                self._config.set(this._section, option, value)
                self._write_config()

            def getboolean(this, option):
                r"""
                TODO
                """
                return self._config.getboolean(this._section, option)

            def __delitem__(this, option):
                r"""
                Remove ``option`` from this section.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: del(c["trac"]["username"])
                    sage: c["trac"]["username"]
                    Traceback (most recent call last):
                    ...
                    KeyError: username

                """
                self._config.remove_option(this._section, option)
                self._write_config()

            def __len__(this):
                r"""
                Return the number of options in this section.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: len(c["trac"])

                """
                return len(self._config.options(this._section))

            def __contains__(this, option):
                r"""
                Return whether this section contains ``options``.

                TESTS::

                    sage: from sage.dev.sagedev.test.config import Config
                    sage: c = Config()
                    sage: "username" in c["trac"]
                    True

                """
                return option in self._config.options(this._section)

        return IndexableForSection(section)

    def __contains__(self, section):
        r"""
        Returns ``True`` if ``section`` is in the configuration.

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: 'trac' in c
            True
            sage: 'nottrac' in c
            False

        """
        return section in self._config.sections()

    def __iter__(self):
        r"""
        Return an iterator over the section names.

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: list(c)
            ['git', 'trac']

        """
        return iter(self._config.sections())

    def __setitem__(self, section, dictionary):
        r"""
        Copy the options from ``dictionary`` to ``section``.

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: c["foo"] = {"foo":"foo"}
            sage: c["foo"]["foo"]
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
        Get the number of sections in the configuration.

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: len(c)
            2

        """
        return len(self._config.sections())

    def __delitem__(self, section):
        r"""
        Remove ``section`` from the configuration

        EXAMPLES::

            sage: from sage.dev.sagedev.test.config import Config
            sage: c = Config()
            sage: del c['git']
            sage: list(c)
            ['trac']

        """
        self._config.remove_section(section)
        self._write_config()
