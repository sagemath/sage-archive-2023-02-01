r"""
Global options

The :class:`AddOptionsToClass` class provides a generic mechanism for
setting and accessing **global** options for parents in one or several
related classes, typically for customizing the representation of their
elements. This class will eventually also support setting options on a
parent by parent basis.

.. SEEALSO::

    For better examples of :class:`AddOptionsToClass` in action see
    :meth:`sage.combinat.partition.Partitions.global_options` and
    :meth:`sage.combinat.tableau.Tableaux.global_options`.

.. _construction_section:

Construction of options classes
-------------------------------

The general setup for creating a set of global options is:

.. code-block:: python

    MyOptions=AddOptionsToClass('option name',
        doc='Nice options',
        first_option=dict(default='default value',
                          description='Changes the functionality of _repr_',
                          values=dict(with_bells='causes _repr_ to print with bells',
                                      with_whistles='causes _repr_ to print with whistles',
                                      ...),
                          aliases=dict(bells='option1', whistles='option2', ...),
        second_option=dict(...),
        third_option=dict(),
        end_doc='end of options documentation'
    )

Each option is specified as a dictionary which describes the possible
values for the option and its documentation. The possible entries in this
dictionary are:

- ``alias`` -- Allows for several option values to do the same thing.

- ``alt_name`` -- An alternative name for this option.

- ``checker`` -- A validation function which returns whether a user
  supplied value is valid or not. This is typically useful for large
  lists of legal values such as `\NN`.

- ``default`` -- Gives the default value for the option.

- ``description`` -- A one line description of the option.

- ``link_to`` -- Links this option to another one in another set of
  global options. This is used for example to allow
  :class:`Partitions` and :class:`Tableaux` to share the same
  ``convention`` option.

- ``setter`` -- A function which is called **after** the value of the
  option is changed.

- ``values`` -- A dictionary assigning each valid value for the option
  to a short description of what it does.

- ``case_sensitive`` -- (Default: ``True``) ``True`` or ``False`` depending on
  whether the values of the option are case sensitive.

For each option, either a complete list of possible values, via ``values``, or a
validation function, via ``checker``, must be given. The values can be quite
arbitrary, including user-defined functions which customize the default
behaviour of the classes such as the output of ``_repr_`` or :func:`latex`. See
:ref:`dispatcher` below, and :meth:`~AddOptionsToClass._dispatcher`, for more
information.

The documentation for the options is automatically constructed by combining the
description of each option with a header and footer which are given by the
following optional, but recommended, arguments:

- ``doc`` -- The top half of the documentation which appears before the
  automatically generated list of options and their possible values.

- ``end_doc`` -- The second half of the documentation which appears
  after the list of options and their values.


The basic structure for defining a :class:`AddOptionsToClass` class is best
illustrated by an example::

    sage: from sage.structure.global_options import AddOptionsToClass
    sage: class MenuClass(object): __name__='Menus'
    sage: AddOptionsToClass(MenuClass, doc='Fancy documentation\n'+'-'*19, end_doc='The END!',
    ...       entree=dict(default='soup',
    ...                   description='The first course of a meal',
    ...                   values=dict(soup='soup of the day', bread='oven baked'),
    ...                   alias=dict(rye='bread')),
    ...       appetizer=dict(alt_name='entree'),
    ...       main=dict(default='pizza', description='Main meal',
    ...                 values=dict(pizza='thick crust', pasta='penne arrabiata'),
    ...                 case_sensitive=False),
    ...       dessert=dict(default='espresso', description='Dessert',
    ...                    values=dict(espresso='life begins again',
    ...                                cake='waist begins again',
    ...                                cream='fluffy, white stuff')),
    ...       tip=dict(default=10, description='Reward for good service',
    ...       checker=lambda tip: tip in range(0,20))
    ...   )
    sage: MenuClass.options()
    options for menu

For more details see :class:`AddOptionsToClass`.

Accessing and setting option values
-----------------------------------

All options and their values, when they are strings, are forced to be lower
case. The values of an options class can be set and accessed by calling the
class or by treating the class as an array.

Continuing the example from :ref:`construction_section`::

    sage: menu()
    Current options for menu
      - dessert: espresso
      - entree:  soup
      - main:    pizza
      - tip:     10
    sage: menu('dessert')
    'espresso'
    sage: menu['dessert']
    'espresso'

Note that, provided there is no ambiguity, options and their values can be
abbreviated::

    sage: menu['d']
    'espresso'
    sage: menu('m','t',des='esp', ent='sou')  # get and set several values at once
    ['pizza', 10]
    sage: menu(t=15); menu['tip']
    15
    sage: menu(e='s', m='Pi'); menu()
    Current options for menu
      - dessert: espresso
      - entree:  soup
      - main:    pizza
      - tip:     15
    sage: menu(m='P')
    Traceback (most recent call last):
    ...
    ValueError: P is not a valid value for main in the options for menu


Setter functions
----------------

Each option of a :class:`AddOptionsToClass` can be equipped with an optional setter
function which is called **after** the value of the option is changed. In the
following example, setting the option 'add' changes the state of the class by
setting an attribute in this class using a :func:`classmethod`. Note that the
options object is inserted after the creation of the class in order to access
the :func:`classmethod` as ``A.setter``::

    sage: from sage.structure.global_options import AddOptionsToClass
    sage: class A(SageObject):
    ...       __name__='A'
    ...       state = 0
    ...       @classmethod
    ...       def setter(cls, option, val):
    ...           cls.state += int(val)
    ...
    sage: AddOptionsToClass(A,
    ...                     add=dict(default=1,
    ...                              checker=lambda v: int(v)>0,
    ...                              description='An option with a setter',
    ...                              setter=A.setter))
    sage: A.options
    sage: a = A(2); a.state
    1
    sage: a.options()
    Current options for A
    - add: 1
    sage: a.options(add=4)
    sage: a.state
    5
    sage: a.options()
    Current options for A
    - add: 4

Another alternative is to construct the options class inside the ``__init__``
method of the class ``A``.

Documentation for options
-------------------------

The documentation for a :class:`AddOptionsToClass` is automatically generated from
the supplied options. For example, the generated documentation for the options
``menu`` defined in :ref:`construction_section` is the following::

    Fancy documentation
    -------------------

    OPTIONS:

    - ``appetizer`` -- alternative name for ``entree``

    - ``dessert`` -- (default: ``espresso``)
      Dessert

      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again

    - ``entree`` -- (default: ``soup``)
      The first course of a meal

      - ``bread`` -- oven baked
      - ``rye``   -- alias for bread
      - ``soup``  -- soup of the day

    - ``main`` -- (default: ``pizza``)
      Main meal

      - ``pasta`` -- penne arrabiata
      - ``pizza`` -- thick crust

    - tip -- (default: 10)
      Reward for good service

    The END!

    See :class:`~sage.structure.global_options.AddOptionsToClass` for more features of these options.

In addition, help on each option, and its list of possible values, can be
obtained by (trying to) set the option equal to '?'::

    sage: menu(des='?')
    - ``dessert`` -- (default: ``espresso``)
      Dessert
    <BLANKLINE>
      - ``cake``     -- waist begins again
      - ``cream``    -- fluffy, white stuff
      - ``espresso`` -- life begins again
    <BLANKLINE>
    Current value: espresso

.. _dispatcher:

Dispatchers
-----------

The whole idea of a :class:`AddOptionsToClass` class is that the options change the
default behaviour of the associated classes. This can be done either by simply
checking what the current value of the relevant option is. Another possibility
is to use the options class as a dispatcher to associated methods. To use the
dispatcher feature of a :class:`AddOptionsToClass` class it is necessary to implement
separate methods for each value of the option where the naming convention for
these methods is that they start with a common prefix and finish with the value
of the option.

If the value of a dispatchable option is set equal to a (user defined) function
then this function is called instead of a class method.

For example, the options ``MyOptions`` can be used to dispatch the ``_repr_``
method of the associated class ``MyClass`` as follows:

.. code-block:: python

    class MyClass(...):
        def _repr_(self):
            return self.options._dispatch(self,'_repr_','first_option')
        def _repr_with_bells(self):
            print 'Bell!'
        def _repr_with_whistles(self):
            print 'Whistles!'
    AddOptionsToClass(MyClass,
        ...
    )

In this example, ``first_option`` is an option of ``MyOptions`` which takes
values ``bells``, ``whistles``, and so on. Note that it is necessary to make
``self``, which is an instance of ``MyClass``, an argument of the dispatcher
because :meth:`~AddOptionsToClass._dispatch()` is a method of :class:`GlobalOptions`
and not a method of ``MyClass``. Apart from ``MyOptions``, as it is a method of
this class, the arguments are the attached class (here ``MyClass``), the prefix
of the method of ``MyClass`` being dispatched, the option of ``MyOptions``
which controls the dispatching. All other arguments are passed through to the
corresponding methods of ``MyClass``. In general, a dispatcher is invoked as::

    self.options._dispatch(self, dispatch_to, option, *args, **kargs)

Usually this will result in the method
``dispatch_to + '_' + MyOptions(options)`` of ``self`` being called with
arguments ``*args`` and ``**kargs`` (if ``dispatch_to[-1] == '_'`` then the
method ``dispatch_to + MyOptions(options)`` is called).

If ``MyOptions(options)`` is itself a function then the dispatcher will call
this function instead. In this way, it is possible to allow the user to
customise the default behaviour of this method. See
:meth:`~AddOptionsToClass._dispatch` for an example of how this can be achieved.

The dispatching capabilities of :class:`AddOptionsToClass` allows options to be
applied automatically without needing to parse different values of the option
(the cost is that there must be a method for each value). The dispatching
capabilities can also be used to make one option control several methods:

.. code-block:: python

    def __le__(self, other):
        return self.options._dispatch(self, '_le_','cmp', other)
    def __ge__(self, other):
        return self.options._dispatch(self, '_ge_','cmp', other)
    def _le_option_a(self, other):
        return ...
    def _ge_option_a(self, other):
        return ...
    def _le_option_b(self, other):
        return ...
    def _ge_option_b(self, other):
        return ...

See :meth:`~AddOptionsToClass._dispatch` for more details.

Doc testing
-----------

All of the options and their effects should be doc-tested. However, in order
not to break other tests, all options should be returned to their default state
at the end of each test. To make this easier, every :class:`AddOptionsToClass` class has
a :meth:`~AddOptionsToClass._reset()` method for doing exactly this.


Tests
-----

TESTS:

As options classes to not know how they are created they cannot be
pickled::

    sage: class MenuClass(object): __name__='Menus'
    sage: AddOptionsToClass(MenuClass, doc='Fancy documentation\n'+'-'*19, end_doc='The END!',
    ...       entree=dict(default='soup',
    ...                   description='The first course of a meal',
    ...                   values=dict(soup='soup of the day', bread='oven baked'),
    ...                   alias=dict(rye='bread')),
    ...       appetizer=dict(alt_name='entree'),
    ...       main=dict(default='pizza', description='Main meal',
    ...                 values=dict(pizza='thick crust', pasta='penne arrabiata'),
    ...                 case_sensitive=False),
    ...       dessert=dict(default='espresso', description='Dessert',
    ...                    values=dict(espresso='life begins again',
    ...                                cake='waist begins again',
    ...                                cream='fluffy, white stuff')),
    ...       tip=dict(default=10, description='Reward for good service',
    ...       checker=lambda tip: tip in range(0,20))
    ...   )
    sage: TestSuite(menu).run()


== Slowest module imports (excluding / including children) ==
exclude/ms include/ms   #parents  module name
     3.593      3.812          1  sage.combinat.posets.hasse_diagram
     4.394     20.471         38  sage.combinat.partition
     4.856      4.935          1  sage.combinat.diagram_algebras
     5.179      6.218         30  sage.combinat.permutation
   208.358    216.462         22  sage.libs.pari.pari_instance
Total time (sum over exclusive time): 1462.549ms
Use sage -startuptime <module_name> to get more details about <module_name>.

AUTHORS:

- Andrew Mathas (2013): initial version
- Andrew Mathas (2016): overhaul making the options attributes, enabling
                        pickling and attaching the options to a class.
"""
#*****************************************************************************
#  Copyright (C) 2013,2016 Andrew Mathas <andrew dot mathas at sydney dot edu dot au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from __builtin__ import object, str
from sage.misc.superseded import deprecated_function_alias
import inspect

class _Option(object):
    r"""
        Each option for an options class is an instance of this class which
        implements the magic that allows the options to the attributes of the
        options class that can be looked up, set and called.

        By way of example, this class implements the following functionality.

        EXAMPLES::

            sage: Partitions.options.display           # indirect doctest
            sage: Partitions.options.display='compact'
            sage: Partitions.options.display('list')

        TESTS::

            sage: TestSuite(Partitions.options.display).run()
    """
    def __init__(self, options, name):
        r"""
        Initialise an option by settings its ``name``, "parent" option class
        ``options`` and doc-string.

        EXAMPLES::

            sage: type(Partitions.options.display)    # indirect doctest
        """
        self._name = name
        self._options = options
        self.__doc__= options._doc[name]
        super(_Option, self).__init__()

    def __repr__(self):
        r"""
        Return a string representation for this collection of options.

        EXAMPLES::

            sage: Partitions.options.display # indirect doctest
        """
        return self._options.__getitem__(self._name)

    def __call__(self, value=None):
        r"""
        Get or set value of the option ``self``.

        EXAMPLES::

            sage: Partitions.options.display() # indirect doctest
            sage: Partitions.options.display('exp') # indirect doctest
        """
        print('Options.__call: {}'.format(value))
        if value is None:
            return self._options[self._name]
        else:
            self._options.__setitem__(self._name, value)

class AddOptionsToClass(object):
    r"""
    The :class:`AddOptionsToClass` class is a generic class for setting and
    accessing global options for ``sage`` objects. It takes as inputs a
    ``name`` for the collection of options and a dictionary of dictionaries
    which specifies the individual options. The allowed/expected keys in the
    dictionary are the following:

    INPUT:

    - ``name`` -- Specifies a name for the options class (required)

    - ``doc`` -- Initial documentation string

    - ``end_doc`` -- Final documentation string

    - ``<options>=dict(...)`` -- Dictionary specifying an option

    The options are specified by keyword arguments with their values
    being a dictionary which describes the option. The
    allowed/expected keys in the dictionary are:

    - ``alias`` -- Defines alias/synonym for option values
    - ``alt_name`` -- Alternative name for an option
    - ``checker`` -- A function for checking whether a particular value for
      the option is valid
    - ``default`` -- The default value of the option
    - ``description`` -- Documentation string
    - ``link_to`` -- Links to an option for this set of options to an
      option in another :class:`AddOptionsToClass`
    - ``setter`` -- A function (class method) which is called whenever this
      option changes
    - ``values`` -- A dictionary of the legal values for this option (this
      automatically defines the corresponding ``checker``). This dictionary
      gives the possible options, as keys, together with a brief description
      of them.
    - ``case_sensitive`` -- (Default: ``True``) ``True`` or ``False`` depending on
      whether the values of the option are case sensitive.

    Options and their values can be abbreviated provided that this
    abbreviation is a prefix of a unique option.

    Calling the options with no arguments results in the list of
    current options being printed.

    EXAMPLES::

        sage: from sage.structure.global_options import AddOptionsToClass
        sage: class Menu(object): __name__='Menus'
        sage: AddOptionsToClass(Menu, doc='Fancy documentation\n'+'-'*19, end_doc='End of Fancy documentation',
        ...       entree=dict(default='soup',
        ...                   description='The first course of a meal',
        ...                   values=dict(soup='soup of the day', bread='oven baked'),
        ...                   alias=dict(rye='bread')),
        ...       appetizer=dict(alt_name='entree'),
        ...       main=dict(default='pizza', description='Main meal',
        ...                 values=dict(pizza='thick crust', pasta='penne arrabiata'),
        ...                 case_sensitive=False),
        ...       dessert=dict(default='espresso', description='Dessert',
        ...                    values=dict(espresso='life begins again',
        ...                                cake='waist begins again',
        ...                                cream='fluffy white stuff')),
        ...       tip=dict(default=10, description='Reward for good service',
        ...                checker=lambda tip: tip in range(0,20))
        ...   )
        sage: Menu.options
        options for menu
        sage: menu(entree='s')         # unambiguous abbreviations are allowed
        sage: menu(t=15);
        sage: (menu['tip'], menu('t'))
        (15, 15)
        sage: menu()
        Current options for menu
          - dessert: espresso
          - entree:  soup
          - main:    pizza
          - tip:     15
        sage: menu._reset(); menu()
        Current options for menu
          - dessert: espresso
          - entree:  soup
          - main:    pizza
          - tip:     10
        sage: menu['tip']=40
        Traceback (most recent call last):
        ...
        ValueError: 40 is not a valid value for tip in the options for menu
        sage: menu(m='p')           # ambiguous abbreviations are not allowed
        Traceback (most recent call last):
        ...
        ValueError: p is not a valid value for main in the options for menu

    The documentation for the options class is automatically generated from the
    information which specifies the options::

        Fancy documentation
        -------------------

        OPTIONS:

        - dessert:  (default: espresso)
          Dessert

          - ``cake``     -- waist begins again
          - ``cream``    -- fluffy white stuff
          - ``espresso`` -- life begins again

        - entree:  (default: soup)
          The first course of a meal

          - ``bread`` -- oven baked
          - ``rye``   -- alias for bread
          - ``soup``  -- soup of the day

        - main:  (default: pizza)
          Main meal

          - ``pasta`` -- penne arrabiata
          - ``pizza`` -- thick crust

        - tip:  (default: 10)
          Reward for good service

        End of Fancy documentation

        See :class:`~sage.structure.global_options.AddOptionsToClass` for more features of these options.

    The possible values for an individual option can be obtained by
    (trying to) set it equal to '?'::

        sage: menu(des='?')
        - ``dessert`` -- (default: ``espresso``)
          Dessert
        <BLANKLINE>
          - ``cake``     -- waist begins again
          - ``cream``    -- fluffy white stuff
          - ``espresso`` -- life begins again
        <BLANKLINE>
        Current value: espresso
    """
    def __init__(self, name='', doc='', end_doc='', **options):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: menu = AddOptionsToClass('menu', doc='Fancy documentation\n'+'-'*19, end_doc='End of Fancy documentation',
            ...       entree=dict(default='soup',
            ...                   description='The first course of a meal',
            ...                   values=dict(soup='soup of the day', bread='oven baked'),
            ...                   alias=dict(rye='bread')),
            ...       appetizer=dict(alt_name='entree'),
            ...       main=dict(default='pizza', description='Main meal',
            ...                 values=dict(pizza='thick crust', pasta='penne arrabiata'),
            ...                 case_sensitive=False),
            ...       dessert=dict(default='espresso', description='Dessert',
            ...                    values=dict(espresso='life begins again',
            ...                                cake='waist begins again',
            ...                                cream='fluffy white stuff')),
            ...       tip=dict(default=10, description='Reward for good service',
            ...                checker=lambda tip: tip in range(0,20))
            ...   )
            sage: specials = AddOptionsToClass('specials menu', doc='More fancy doc...',
            ...       entree=dict(link_to=(menu, 'entree')),
            ...       main_specials=dict(default='salmon', description='main course specials',
            ...                     values=dict(salmon='a fish', crab='Sebastian'))
            ...   )
            sage: specials['entree'] = 'rye'
            sage: menu['entree']
            'bread'

            sage: alias_test = AddOptionsToClass( name='alias_test',
            ...         doc="Test aliases with case sensitivity",
            ...         test_opt=dict(default="Upper",
            ...         description='Starts with an uppercase',
            ...         values=dict(Upper="Starts with uppercase",
            ...                     lower="only lowercase"),
            ...         case_sensitive=False,
            ...         alias=dict(UpperAlias="Upper", lower_alias="lower")) )
            sage: alias_test['test_opt'] = 'Lower_Alias'
            sage: alias_test['test_opt']
            'lower'
            sage: alias_test['test_opt'] = 'upperalias'
            sage: alias_test['test_opt']
            'Upper'
        """
        # initialise the various dictionaries used by AddOptionsToClass
        self._alias={}            # alias for the values of some options
        self._alt_names={}        # alternative names for some options
        self._checker={}          # validity checkers for each option
        self.__default_value={}   # the default options
        self._doc={}              # the linked options force us to keep a dictionary of doc strings
        self._linked_value={}     # linked to other global options as (link, linked_option)
        self._setter={}           # a dictionary of the list of setters
        self._value={}            # the current options
        self._legal_values={}     # a dictionary of lists of the legal values for each option
        self._display_values={}   # a dictionary of the output of the values
        self._case_sensitive = {} # a dictionary of booleans indicating to check case sensitivity
        for option in options:
            self._add_option(option, options[option])

        # Finally, we build the doc string for the options
        # First we strip common white space off the front of doc and end_doc
        if len(doc)>0:
            lines=doc.splitlines()
            m=min(len(line)-len(line.lstrip()) for line in lines if len(line)>0)
            self._doc_start='\n'.join(line[m:] for line in lines)

        if len(end_doc)>0:
            lines=end_doc.splitlines()
            m=min(len(line)-len(line.lstrip()) for line in lines if len(line)>0)
            self._doc_end='\n'.join(line[m:] for line in lines)

        super(AddOptionsToClass, self).__init__()

    __name__ = 'Options class'

    def __repr__(self):
        r"""
        Return a string representation for this collection of options.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...          food=dict(default='apple', values=dict(apple='a nice fruit',pear='fruit')),
            ...          drink=dict(default='water', values=dict(water='wet',milk='white')))
            sage: FoodOptions
            options for daily meal
        """
        options=self._value.keys()+self._linked_value.keys()
        for x in self._alt_names:
            options.remove(x)
        if options == []:
            return  'Current options for {}'.format(self._name)

        options.sort()
        width=1+max(len(option) for option in options)
        return  'Current options for {}\n{}'.format(self._name,
                    '\n'.join('  - {:{}} {}'.format(option+':',width,self[option]) for option in options)
                )

    def __call__(self, *get_value, **set_value):
        r"""
        Get or set value of the option ``option``.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',coffee='a lifestyle')),
            ...         beverage=dict(alt_name='drink'))
            sage: FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  apple
            sage: FoodOptions('food')
            'apple'
            sage: FoodOptions(food="pair"); FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  pair
            sage: FoodOptions('beverage')
            'water'
            sage: FoodOptions(food="apple", drink="coffee"); FoodOptions()
            Current options for daily meal
              - drink: coffee
              - food:  apple
        """
        if get_value==() and set_value=={}:
            print self
            return

        if get_value!=():
        # use __getitem__ to return these options
            if len(get_value)==1:
                return self.__getitem__(get_value[0])
            else:
                return [self.__getitem__(option) for option in get_value]

        # use __setitem__ to set these options
        if set_value!=[]:
            for option in set_value:
                self.__setitem__(option, set_value[option])

    def __getitem__(self, option):
        r"""
        Return the current value of the option ``option``.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',coffee='a lifestyle')),
            ...         beverage=dict(alt_name='drink'))
            sage: FoodOptions['drink']
            'water'
            sage: FoodOptions['d']
            'water'
        """
        option=self._match_option(option)
        if option in self._linked_value:
            link,linked_opt=self._linked_value[option]
            return link[linked_opt]
        elif option in self._value:
            if option in self._display_values:
                return self._display_values[option][self._value[option]]
            return self._value[option]

    def __setitem__(self, option, value):
        r"""
        The ``__setitem__`` method is used to change the current values of the
        options. It also checks that the supplied options are valid and changes
        any alias to its generic value.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',coffee='a lifestyle')))
            sage: FoodOptions['drink']='coffee'; FoodOptions()
            Current options for daily meal
              - drink: coffee
              - food:  apple
            sage: FoodOptions(drink='w'); FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  apple
            sage: FoodOptions(drink='?')
            - ``drink`` -- (default: ``water``)
            <BLANKLINE>
              - ``coffee`` -- a lifestyle
              - ``water``  -- a drink
            <BLANKLINE>
            Current value: water
        """
        option=self._match_option(option)
        if not callable(value):
            value=self._match_value(option, value)

        if value=='?':  # return help
            print('%s\nCurrent value: %s' % (self._doc[option], self[option]))
            return      # we do not want to call the setter below

        elif option in self._linked_value:
            link, linked_opt=self._linked_value[option]
            link[linked_opt]=value

        else:
            self._value[option]=value

        if option in self._setter:
            # if a setter function exists then call it with the associated
            # class, option and value
            self._setter[option](option, value)

    def __docstring__ (self):
        r''' Return the docstring for the options class ``self``.

            EXAMPLES::

                sage: Partitions.options?  # not tested (too long)
        '''
        return '{start}\n\nOPTIONS:\n\n{options}\n\n\n{end}\n\n{g_opts}'.format(
                   start=self._doc_start, end=self._doc_end,
                   options='\n'.join(self._doc[opt] for opt in sorted(self._doc)),
                   g_opts='See :class:`~sage.structure.global_options.AddOptionsToClass` for more features of these options.'
        )

    def __getattribute__(self, name):
        r"""
        Return the attribute ``name`` for ``self.

        If ``name`` is `__doc__` then we return ``self._docstring``. This allows
        us to dynamically alter the doc string when new options are added.

        EXAMPLES::

            sage: Partitions.options.display # indirect docttest
        """
        if name == '__doc__':
            return object.__getattribute__(self, '__docstring__')()
        else:
            return object.__getattribute__(self, name)

    def __setattr__(self, name, value=None):
        r"""
            Set the attribute ``name`` of the option class self equal to ``value,  
            if the attribute ``name`` exists.

            As the attributes of an option class are the actual options we need
            to be able to "trap" invalid options in a sensible way. We do this
            by sending any "non-standard" to :meth:`__setitem__` for processing.

            EXAMPLES::

                sage: Partitions.options.display='exp'  # indirect doc-test
                sage: Partitions.options.dispplay='exp' # indirect doc-test
        """
        # Underscore, and "special", attributes are set using object.__setattr__
        # Anything else is assume to be an option and directed to __setitem__.
        if name[0] == '_' or name in ['reset', 'dispatch', 'default_value']:
            object.__setattr__(self, name, value)
        else: # redirect to __setitem
            self.__setitem__(name, value)

    def __setstate__(self, state):
        r"""
        This is a custom :meth:`__setstate__` method for unpickling instances of
        the :class:`AddOptionsToClass` class.

        The :meth:`__getstate__` method returns a dictionary with an
        `options_class` key which identifies the "parent" class for the options.
        This is then used to unpickle the options class.

        EXAMPLES::

            sage: Partitions.global_options()
            Current options for Partitions
              - convention:        English
              - diagram_str:       *
              - display:           list
              - latex:             young_diagram
              - latex_diagram_str: \ast
            sage: Partitions.options.convention="French"
            sage: loads(dumps(Partitions.global_options))()  # indirect doctest
            Current options for Partitions
              - convention:        French
              - diagram_str:       *
              - display:           list
              - latex:             young_diagram
              - latex_diagram_str: \ast
        """
        # copy all settings across from unpickle to `self`.
        unpickle=state['options_class'].options
        for setting in unpickle.__dict__.keys():
            self.__dict__[setting] = unpickle.__dict__[setting]
        self._reset()       # reset the options in `self` to their defaults
        state.pop('options_class')
        for opt in state: # apply the options store in state
            self[opt]=state[opt]
        self._options_class.options=self

    def __getstate__(self):
        r"""
        Returns a dictionary that can be used to pickle an instance of a
        :class:`AddOptionsToClass`  class.

        Typically instances of :class:`AddOptionsToClass` are complicated to create
        so they do no pickle. If the options are associated to  "parent" class
        then it is possible to create the default version of the class and then
        add the non-default settings on top of this. This method returns a
        dictionary of the non-default options that can then be used to unpickle
        an instance of the option using :meth:`__setstate__`.

        EXAMPLES::

            sage: Partitions.options._reset()
            sage: Partitions.options.__getstate__()
                {'convention': 'English',
                 'options_class': <class 'sage.combinat.partition.Partitions'>}
        """
        try: 
            options_class=globals()[self._name]
            if inspect.isclass(option_class) and hasattr(self._options_class, 'options'):
                pickle={'options_class': self._options_class}
                for opt in self._value.keys():
                    if opt not in self._alt_names and self[opt]!=self.__default_value[opt]:
                        pickle[opt]=self[opt]
                for opt in self._linked_value:
                    link, linked_opt=self._linked_value[opt]
                    if opt not in self._alt_names and link[opt]!=link.__default_value[opt]:
                        pickle[opt]=self[opt]

                try:
                    return pickle
                except PicklingError:
                    raise PicklingError('one or more of the options for %s cannot  be pickled' % self._options_class)
        except KeyError:
            pass

        # if self._options_class is not a class then we have no way to
        # reconstruct the global options
        raise PicklingError('%s cannot be pickled because it is not associated with a class' % self)

    def __eq__(self, other):
        r"""
        Two options classes are equal if they return the same :meth:`__getstate__.

        EXAMPLES::

            sage: Partitions.options == PartitionsGreatestLE.options # indirect doctest
            True
            sage: Partitions.options == Tableaux.options
            False
        """
        return self.__getstate__() == other.__getstate__()

    def _add_option(self, option, specifications):
        r"""
        Add an option.

        INPUT:

        - ``option`` -- a string
        - ``specifications`` -- a dictionary

        .. SEEALSO::

            :class:`AddOptionsToClass` for a description of the required
            ``specifications``.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions = AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',coffee='a lifestyle')),
            ...         beverage=dict(alt_name='drink')) # indirect doctest
            sage: FoodOptions()
            Current options for daily meal
              - drink: water
              - food:  apple
        """
        doc={}  # will be used to build the doc string
        option = option.lower()
        self._legal_values[option] = []
        self._case_sensitive[option] = True    # ``True`` by default
        for spec in sorted(specifications):   # NB: options processed alphabetically!
            if spec=='alias':
                self._alias[option]=specifications[spec]
                self._legal_values[option]+=specifications[spec].keys()
                for opt in specifications[spec]:
                    doc[opt] = 'alias for ``%s``'%specifications[spec][opt]
            elif spec == 'alt_name':
                self._alt_names[option] = specifications[spec]
                self._linked_value[option] = (self, specifications[spec])
                doc = '- ``%s`` -- alternative name for ``%s``'%(option, specifications[spec].lower())
            elif spec=='checker':
                if not callable(specifications[spec]):
                    raise ValueError('the checker for %s must be callable'%option)
                self._checker[option]=specifications[spec]
            elif spec=='default':
                self.__default_value[option]=specifications[spec]
            elif spec=='link_to':
                if (isinstance(specifications[spec], tuple) and len(specifications[spec]) == 2 and
                        isinstance(specifications[spec][0], AddOptionsToClass)):
                    link, linked_opt = specifications['link_to']  # for sanity
                    if linked_opt in link._value:
                        self._linked_value[option] = specifications['link_to']
                        link, linked_opt = specifications['link_to']  # for sanity
                        doc = link._doc[linked_opt]
                    elif linked_opt in link._linked_value:
                        self._linked_value[option] = link._linked_value[linked_opt]
                        doc = link._doc[linked_opt]
                    else:
                        raise ValueError("could not find link to {1} in {0}".format(*specifications[spec]))
                else:
                    raise ValueError("linked options must be specified as a string: 'linked_option' or a tuple: (link,linked_option)")
            elif spec=='setter':
                if callable(specifications[spec]):
                    self._setter[option]=specifications[spec]
                else:
                    raise ValueError('the setter for %s must be a function' % option)
            elif spec=='values':
                for val in specifications[spec]:
                    doc[val] = specifications[spec][val]
                doc.update(specifications[spec])
                if self._case_sensitive[option]:
                    self._legal_values[option] += [val for val in specifications[spec].keys()]
                    self._display_values[option] = {val:val for val in specifications[spec].keys()}
                else:
                    self._legal_values[option] += [val.lower() for val in specifications[spec].keys()]
                    self._display_values[option] = {val.lower():val for val in specifications[spec].keys()}
            elif spec == 'case_sensitive':
                if not specifications[spec]:
                    for opt in self._legal_values:
                        self._display_values[option] = {val.lower():val for val in self._legal_values[option]}
                        self._legal_values[option] = [val.lower() for val in self._legal_values[option]]
                    if option in self._alias:
                        self._alias[option] = {k.lower():v.lower() for k,v in self._alias[option].iteritems()}
                self._case_sensitive[option] = bool(specifications[spec])
            elif spec!='description':
                raise ValueError('Initialization error in Global options for %s: %s not recognized!'%(self._name, spec))

        # now build the doc string for this option
        if doc == {} and not 'description' in specifications:
            raise ValueError('no documentation specified for %s in the %s' % (option, self))

        # first a necessary hack to initialise the option in self._doc because __setitem__ calls _match_option
        self._doc[option]=''   
        if option in self._linked_value:
            self._doc[option]=doc
        else:
            width = max(len(v) for v in doc.keys()) + 4 if doc!={} else 4
            if len(doc) > 0:
                self._doc[option.lower()]='- ``{}`` -- (default: ``{}``)\n{}\n{}\n'.format(
                    option, self._default_value(option),
                    '  %s\n'%specifications['description'] if 'description' in specifications else '',
                    '\n'.join('  - {:{}} -- {}'.format('``'+val+'``',width,doc[val]) for val in sorted(doc)))
            else:
                self._doc[option.lower()]='- ``{}`` -- (default: ``{}``)\n{}'.format(
                    option, self._default_value(option),
                    '  %s\n'%specifications['description'] if 'description' in specifications else '')

        # sanity check for non-linked options
        if not option in self._linked_value:
            if 'default' not in specifications:
                raise ValueError('a default value for %s must be given' % option)

            if not (option in self._checker or option in self._legal_values):
                raise ValueError('a value checker or a list of valid values for %s must be given' % option)

            # finally, set, check and process the default value using  __setitem__
            self[option]=self.__default_value[option]
            self.__default_value[option]=self._value[option]  # in case the default is an alias

        # Build getters and setters for this option. As we have overridden __setattr__ we 
        # call object.__setattr_ directly
        object.__setattr__(self, option, _Option(self, option))

    def _match_option(self, option):
        r"""
        Check whether ``option`` is the name of an option of ``self``, or a
        prefix thereof.

        INPUT:

        - ``option`` -- a string

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',coffee='a lifestyle')))
            sage: FoodOptions('food') # indirect doctest
            'apple'
            sage: FoodOptions('f')
            'apple'
        """
        # the keys of self._doc is a list of the options, both normal and linked
        option = option.lower()

        if option in self._doc: return option

        # as it is not an option try and match it with a prefix to an option
        matches=[opt for opt in self._doc if opt.startswith(option)]
        if len(matches)>0 and all(m.startswith(matches[0]) for m in matches):
            return matches[0]
        elif len(matches)>1:
            raise ValueError('%s is an ambiguous option for %s'%(option, self._name))

        # if we are still here this is not a good option!
        raise ValueError('%s is not an option for %s' % (option, self._name))

    def _match_value(self, option, value):
        r"""
        Check whether ``value`` is a valid value for ``option``.

        INPUT:

        - ``option`` -- a string: the name of an option, or prefix thereof
        - ``value``  -- a value or ``'?'``

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',wine='a lifestyle')))
            sage: FoodOptions(f='a') # indirect doctest
            sage: FoodOptions('f')
            'apple'
            sage: FoodOptions(d='wi'); FoodOptions('f')
            'apple'
            sage: FoodOptions(d='w')
            Traceback (most recent call last):
            ...
            ValueError: w is not a valid value for drink in the options for daily meal
        """
        if value == "?": return value   # help on this value

        if option in self._linked_value:
            link, linked_opt = self._linked_value[option]
            return link._match_value(linked_opt, value)

        orig_value = value

        # convert to proper case if it is a string
        if isinstance(value, str) and not self._case_sensitive[option]:
            value = value.lower()

        if option in self._legal_values:
            if value in self._legal_values[option]:
                if option in self._alias and value in self._alias[option]:
                    return self._alias[option][value]
                return value

            # as it is not a value try and match it with a prefix of a value
            matches=[val for val in self._legal_values[option] if val.startswith(value)]
            if len(matches)>0 and all(m.startswith(matches[0]) for m in matches):
                val=matches[0]
                if option in self._alias and val in self._alias[option]:
                    return self._alias[option][val]
                return val

        if option in self._checker and self._checker[option](value):
            return value

        # if we are still here this is not a good value!

        # replace any value alias with its "real" value
        if option in self._alias and value in self._alias[option]:
            orig_value = self._alias[option][value]
        raise ValueError('%s is not a valid value for %s in the %s'%(orig_value, option, self))

    def _default_value(self, option):
        r"""
        Return the default value of the option.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: FoodOptions=AddOptionsToClass('daily meal',
            ...         food=dict(default='apple', values=dict(apple='a fruit',pair='of what?')),
            ...         drink=dict(default='water', values=dict(water='a drink',wine='a lifestyle')))
            sage: FoodOptions._default_value('food')
            'apple'
        """
        option=self._match_option(option)
        if option in self.__default_value:
            return self.__default_value[option]
        else:
            link, linked_opt=self._linked_value[option]
            return link._default_value(linked_opt)

    default_value=deprecated_function_alias(18555, _default_value)

    def _dispatch(self, obj, dispatch_to, option, *args, **kargs):
        r"""
        .. TODO:: title

        The *dispatchable* options are options which dispatch related methods of
        the corresponding class - or user defined methods which are passed to
        :class:`AddOptionsToClass`. The format for specifying a dispatchable option
        is to include ``dispatch_to = <option name>`` in the specifications for
        the options and then to add the options to the (element) class. Each
        option is then assumed to be a method of the element class with a name
        of the form ``<option name> + '_' + <current vale for this option'``.
        These options are called by the element class via::

            return self.options._dispatch(self, dispatch_to, option, *args, **kargs)

        Note that the argument ``self`` is necessary here because the
        dispatcher is a method of the options class and not of
        ``self``. The value of the option can also be set to a
        user-defined function, with arguments ``self`` and ``option``;
        in this case the user's function is called instead.

        EXAMPLES:

        Here is a contrived example::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: DelimitedListOptions=AddOptionsToClass('list delimiters',
            ...             delim=dict(default='b', values={'b':'brackets', 'p':'parentheses'}))
            sage: class DelimitedList(CombinatorialObject):
            ...      options=DelimitedListOptions
            ...      def _repr_b(self): return '[%s]' % ','.join('%s'%i for i in self._list)
            ...      def _repr_p(self): return '(%s)' % ','.join('%s'%i for i in self._list)
            ...      def _repr_(self): return self.options._dispatch(self, '_repr_','delim')
            sage: dlist=DelimitedList([1,2,3]); dlist
            [1,2,3]
            sage: dlist.options(delim='p'); dlist
            (1,2,3)
            sage: dlist.options(delim=lambda self: '<%s>' % ','.join('%s'%i for i in self._list)); dlist
            <1,2,3>
        """
        if callable(self._value[option]):
            try:
                return self._value[option](obj, *args, **kargs)
            except TypeError:
                raise ValueError('the user defined dispatcher function failed!')
        else:
            if dispatch_to[-1]=='_': dispatch_to=dispatch_to[:-1]
            dispatch=getattr(obj, dispatch_to+'_'+self._value[option])
            return dispatch(*args, **kargs)

        raise ValueError('%s is not a dispatchable option!' % option)

    dispatch=deprecated_function_alias(18555, _dispatch)

    def _reset(self, option=None):
        r"""
        Reset options to their default value.

        INPUT:

        - ``option`` -- (Default: ``None``) The name of an option as a string
          or ``None``. If ``option`` is specified only this option is reset to
          its default value; otherwise all options are reset.

        EXAMPLES::

            sage: from sage.structure.global_options import AddOptionsToClass
            sage: class Meal(object): __name__='daily meal'
            sage: AddOptionsToClass('daily meal',
            ...      food=dict(default='bread', values=dict(bread='rye bread', salmon='a fish')),
            ...      drink=dict(default='water',values=dict(water='essential for life',wine='essential')))
            sage: Meal.options(food='salmon'); opts()
            Current options for daily meal
              - drink: water
              - food:  salmon
            sage: Meal.options._reset('drink'); Meal.options()
            Current options for daily meal
              - drink: water
              - food:  salmon
            sage: Meal.options._reset(); Meal.options()
            Current options for daily meal
              - drink: water
              - food:  bread
        """
        if option is None:
            for option in self.__default_value:
                self._value[option] = self.__default_value[option]
                if not self._case_sensitive[option] and isinstance(self._value[option],str):
                    self._value[option] = self._value[option].lower()
            for option in self._linked_value:
                link, linked_opt=self._linked_value[option]
                link._reset(linked_opt)
        else:
            option=self._match_option(option)
            if option in self.__default_value:
                self._value[option] = self.__default_value[option]
                if not self._case_sensitive[option] and isinstance(self._value[option],str):
                    self._value[option] = self._value[option].lower()
            elif option in self._linked_value:
                link, linked_opt=self._linked_value[option]
                link._reset(linked_opt)

    reset=deprecated_function_alias(18555, _reset)

GlobalOptions = deprecated_function_alias(18555, AddOptionsToClass)
