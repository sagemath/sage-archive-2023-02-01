# -*- coding: utf-8 -*-
r"""
Factory for Character-Based Art
"""
#*******************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
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
#*******************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.all import ZZ


class CharacterArtFactory(SageObject):

    def __init__(self,
                 art_type, string_type, magic_method_name,
                 parenthesis, square_bracet, curly_brace):
        """
        Abstract base class for character art factory

        This class is the common implementation behind
        :func:`~sage.typeset.ascii_art.ascii_art` and
        :func:`sage.typeset.unicode_art.unicode_art`.

        INPUT:

        - ``art_type` -- type of the character art.

        - ``string_type`` -- type of strings (the lines in the
          character art).

        - ``magic_method_name`` -- name of the Sage magic method.

        - ``parenthesis`` -- left/right pair of two multi-line
          symbols. The parenthesis, a.k.a. round brackets.

        - ``square_bracket`` -- left/right pair of two multi-line
          symbols. The square_brackets.

        - ``curly_brace`` -- left/right pair of two multi-line
          symbols. The curly braces.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: type(factory)
            <class 'sage.typeset.character_art_factory.CharacterArtFactory'>
        """
        self.art_type = art_type
        assert string_type in [str, unicode]
        self.string_type = string_type
        assert magic_method_name in ['_ascii_art_', '_unicode_art_']
        self.magic_method_name = magic_method_name
        self.left_parenthesis, self.right_parenthesis = parenthesis
        self.left_square_bracket, self.right_square_bracket = square_bracet
        self.left_curly_brace, self.right_curly_brace = curly_brace
        
    def build(self, obj):
        r"""
        Construct a character art reprensentation

        INPUT:

        - ``obj`` -- anything. The object whose ascii art representation
          we want.

        OUTPUT:
        
        Character art object.

        EXAMPLES::
    
            sage: ascii_art(integral(exp(x+x^2)/(x+1), x))
                /
               |
               |   2
               |  x  + x
               | e
               | ------- dx
               |  x + 1
               |
              /
    
        TESTS::
    
            sage: n = var('n')
            sage: ascii_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
             /        __________    \
            -\2*x + \/ -4*x + 1  - 1/
            --------------------------
                       __________
                 2*x*\/ -4*x + 1
            sage: ascii_art(list(DyckWords(3)))
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: ascii_art(1)
            1
        """
        if isinstance(obj, self.art_type):
            return obj
        elif isinstance(obj, (tuple, list, dict, set)):
            if obj.__class__ is tuple:
                return self.build_tuple(obj)
            elif obj.__class__ is dict:
                return self.build_dict(obj)
            elif obj.__class__ is list:
                return self.build_list(obj)
            else:
                return self.build_set(obj)
        elif isinstance(obj, SageObject):
            return self.build_from_magic_method(obj)
        else:
            return self.build_from_string(obj)
    
    def build_empty(self):
        """
        Return the empty character art object

        OUTPUT:

        Character art instance.
        """
        return self.art_type.empty()

    def build_from_magic_method(self, obj):
        """
        Return the character art object created by the object's magic method

        OUTPUT:

        Character art instance.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: out = factory.build_from_magic_method(identity_matrix(2));  out
            [1 0]
            [0 1]
            sage: type(out)
            <class 'sage.typeset.ascii_art.AsciiArt'>
        """
        magic_method = getattr(obj, self.magic_method_name)
        return magic_method()

    def build_from_string(self, obj):
        """
        Return the character art object created from splitting the object's string representation

        OUTPUT:

        Character art instance.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: out = factory.build_from_string('a\nbb\nccc')
            sage: out + out + out
            a  a  a
            bb bb bb
            ccccccccc
            sage: type(out)
            <class 'sage.typeset.ascii_art.AsciiArt'>
        """
        lines = self.string_type(obj).splitlines()
        return self.art_type(lines)
    
    def build_comma_sequence(self, iterable, func=None):
        r"""
        Build up a comma-separated sequence

        Auxiliary function for ``build_X`` where ``X`` is ``dict``,
        ``set``, ``list``, or ``tuple``.
    
        EXAMPLES::
    
            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: out = factory.build_comma_sequence(list(DyckWords(3)));  out
                                              /\  
                       /\    /\      /\/\    /  \ 
            /\/\/\, /\/  \, /  \/\, /    \, /    \
            sage: type(out)
            <class 'sage.typeset.ascii_art.AsciiArt'>
        """
        if func is None:
            func = lambda elem, _: self.build(elem)
        repr_elems = self.build_empty()
        not_first = False
        for elem in iterable:
            if not_first:
                #print repr_elems.get_baseline()
                if repr_elems._baseline is not None:
                    repr_elems += self.art_type(
                        ["" ** ZZ(repr_elems._h - 1 - repr_elems._baseline) ] +
                        [", "],
                        baseline=repr_elems._baseline)
                else:
                    repr_elems += self.art_type(
                        (["" ** ZZ(repr_elems._h - 1) ] if repr_elems._h > 1 else [])+
                        [", "],
                        baseline=repr_elems._baseline)
                repr_elems._breakpoints.append(repr_elems._l - 1)
            repr_elems += func(elem, iterable)
            not_first = True
        return repr_elems
    
    def build_set(self, set):
        r"""
        Return an character art output of a set.
    
        TESTS::
    
            sage: ascii_art(set(DyckWords(3)))
            {                                   /\   }
            {  /\      /\/\              /\    /  \  }
            { /  \/\, /    \, /\/\/\, /\/  \, /    \ }
        """
        repr_elems = self.build_comma_sequence(set)
        return self.build_container(
            repr_elems, self.left_curly_brace, self.right_curly_brace)
    
    def build_dict(self, dict):
        r"""
        Return an character art output of a dictionnary.
    
        TESTS::
    
            sage: ascii_art({i:dw for i,dw in enumerate(DyckWords(3))})
            {                                             /\   }
            {                /\      /\        /\/\      /  \  }
            { 0:/\/\/\, 1:/\/  \, 2:/  \/\, 3:/    \, 4:/    \ }
        """
        colon = self.art_type([self.string_type(':')], baseline=0)
        def func(k, dict):
            return self.build(k) + colon + self.build(dict[k])
        repr_elems = self.build_comma_sequence(dict, func)
        return self.build_container(
            repr_elems, self.left_curly_brace, self.right_curly_brace)
    
    def build_container(self, content, left_border, right_border):
        r"""
        Return character art for a container
    
        INPUT:

        - ``content`` --
          :class:`~sage.typeset.character_art.CharacterArt`. The
          content of the container, usually comma-separated entries.

        - ``left_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`. The left
          border of the container.

        - ``right_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`. The right
          border of the container.

        TESTS::
    
            sage: l = ascii_art(list(DyckWords(3)));  l
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: l.get_breakpoints()
            [2, 8, 9, 10, 16, 17, 18, 24, 25, 26, 32, 33, 34, 40]
        """
        # new_mat = []
        # for line in repr_elems:
        #    new_mat.append(l_border + line + " "**ZZ(len(repr_elems) - len(line)) + r_border)
        w = content.width()
        h = content.height()
        left_border = left_border.character_art(h)
        right_border = right_border.character_art(h)
        lines = []
        pad = self.string_type(' ')
        for left, line, right in zip(left_border, content, right_border):
            lines.append(left + pad + line.ljust(w) + pad + right)
        shift = len(left_border) + len(pad)
        basepoints = [bp + shift for bp in content.get_breakpoints()] + [w+shift]
        return self.art_type(lines, basepoints, baseline=0, atomic=False)
    
    def build_list(self, list):
        r"""
        Return an character art output of a list.
    
        TESTS::
    
            sage: l = ascii_art(list(DyckWords(3)));  l
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: l.get_breakpoints()
            [2, 8, 9, 10, 16, 17, 18, 24, 25, 26, 32, 33, 34, 40]
        """
        repr_elems = self.build_comma_sequence(list)
        return self.build_container(
            repr_elems, self.left_square_bracket, self.right_square_bracket)
    
    def build_tuple(self, tuple):
        r"""
        Return an character art output of a tuple.
    
        TESTS::
    
            sage: ascii_art(tuple(DyckWords(3)))
            (                                   /\   )
            (            /\    /\      /\/\    /  \  )
            ( /\/\/\, /\/  \, /  \/\, /    \, /    \ )
        """
        repr_elems = self.build_comma_sequence(tuple)
        return self.build_container(
            repr_elems, self.left_parenthesis, self.right_parenthesis)
    
    def concatenate(self, iterable, separator, empty):
        """
        Concatenate multiple character art instances

        INPUT:

        - ``iterable`` -- iterable of character art.

        - ``separable`` -- character art. The separator in-between the
          iterable.

        - ``empty`` -- the empty character art.

        EXAMPLES::

            sage: i2 = identity_matrix(ZZ, 2)
            sage: ascii_art(i2, i2, i2, sep=ascii_art(1/x))
                 1     1
            [1 0]-[1 0]-[1 0]
            [0 1]x[0 1]x[0 1]
        """
        result = empty
        for obj in iterable:
            if result is not empty:
                result += separator
            result += obj
        return result
