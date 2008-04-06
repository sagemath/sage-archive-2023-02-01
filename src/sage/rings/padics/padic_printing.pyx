"""
This file contains code for printing p-adic elements.

It has been moved here to prevent code duplication and make finding the relevant code easier.
"""

include "../../ext/stdsage.pxi"
include "../../ext/gmp.pxi"

import sys

from sage.rings.integer cimport Integer

cdef enum print_modes:
    terse
    series
    val_unit
    digits
    bars


def pAdicPrinter(ring, mode = None, pos = None, pname = None, unram_name = None, var_name = None, max_ram_terms = None, max_unram_terms = None, max_terse_terms = None, sep = None, alphabet = None):
    return pAdicPrinter_class(ring, mode, pos, pname, unram_name, var_name, max_ram_terms, max_unram_terms, max_terse_terms, sep, alphabet)

class pAdicPrinterDefaults(SageObject):

    def __init__(self, mode = 'series', pos = True, max_ram_terms = -1, max_unram_terms = -1, max_terse_terms = -1, sep = "|", alphabet = None):
        self._mode = mode
        self._pos = bool(pos)
        if not -1 <= max_ram_terms <= sys.maxint:
            raise ValueError, "max_ram_terms must be positive and fit in a long"
        self._max_ram_terms = int(max_ram_terms)
        if not -1 <= max_unram_terms <= sys.maxint:
            raise ValueError, "max_unram_terms must be positive and fit in a long"
        self._max_unram_terms = int(max_unram_terms)
        if not -1 <= max_terse_terms <= sys.maxint:
            raise ValueError, "max_terse_terms must be positive and fit in a long"
        self._max_terse_terms = int(max_terse_terms)
        self._sep = sep
        if alphabet is None:
            self._alphabet = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
        else:
            self._alphabet = alphabet

    def mode(self, mode=None):
        """
        Set the default printing mode.

        EXAMPLES:
            sage: padic_printing.mode('terse')
            sage: padic_printing.mode()
            'terse'
            sage: Qp(7)(100)
            100 + O(7^20)
            sage: padic_printing.mode('series')
            sage: Qp(11)(100)
            1 + 9*11 + O(11^20)
            sage: padic_printing.mode('val-unit')
            sage: Qp(13)(130)
            13 * 10 + O(13^21)
            sage: padic_printing.mode('digits')
            sage: repr(Qp(17)(100))
            '...5F'
            sage: repr(Qp(17)(1000))
            '...37E'
            sage: padic_printing.mode('bars')
            sage: repr(Qp(19)(1000))
            '...2|14|12'

            sage: padic_printing.mode('series')
        """
        if mode is None:
            return self._mode
        else:
            if mode in ['val-unit','series','terse','digits','bars']:
                self._mode = mode
            else:
                raise ValueError, "invalid printing mode"

    def allow_negatives(self, neg = None):
        """
        Controls whether or not to display a balanced representation.

        EXAMPLES:
            sage: padic_printing.allow_negatives(True)
            sage: padic_printing.allow_negatives()
            True
            sage: Qp(29)(-1)
            -1 + O(29^20)
            sage: Qp(29)(-1000)
            -14 - 5*29 - 29^2 + O(29^20)
            sage: padic_printing.allow_negatives(False)
        """
        if neg is None:
            return not self._pos
        else:
            self._pos = not neg

    def max_series_terms(self, max = None):
        """
        EXAMPLES:
            sage: padic_printing.max_series_terms(2)
            sage: padic_printing.max_series_terms()
            2
            sage: Qp(31)(1000)
            8 + 31 + ... + O(31^20)
            sage: padic_printing.max_series_terms(-1)
            sage: Qp(37)(100000)
            26 + 37 + 36*37^2 + 37^3 + O(37^20)
        """
        if max is None:
            return self._max_ram_terms
        else:
            self._max_ram_terms = int(max)

    def max_unram_terms(self, max = None):
        """
        EXAMPLES:
            sage: padic_printing.max_unram_terms(2)
            sage: padic_printing.max_unram_terms()
            2
            sage: Zq(5^6, 5, names='a')([1,2,3,-1])^17
            (2*a + 3 + ...) + (a^2 + a + ...)*5 + (a + 2 + ...)*5^2 + (a + 2 + ...)*5^3 + (4*a + 4 + ...)*5^4 + O(5^5)

            sage: padic_printing.max_unram_terms(-1)
        """
        if max is None:
            return self._max_unram_terms
        else:
            self._max_unram_terms = int(max)

    def max_poly_terms(self, max = None):
        """
        EXAMPLES:
            sage: padic_printing.max_poly_terms(3)
            sage: padic_printing.max_poly_terms()
            3
            sage: padic_printing.mode('terse')
            sage: Zq(7^5, 5, names='a')([2,3,4])^8
            2570 + 15808*a + 9018*a^2 + ... + O(7^5)

            sage: padic_printing.max_poly_terms(-1)
            sage: padic_printing.mode('series')
        """
        if max is None:
            return self._max_terse_terms
        else:
            self._max_terse_terms = int(max)

    def sep(self, sep = None):
        """
        EXAMPLES:
            sage: padic_printing.sep('][')
            sage: padic_printing.sep()
            ']['
            sage: padic_printing.mode('bars')
            sage: repr(Qp(61)(-1))
            '...60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60][60'

            sage: padic_printing.sep('|')
            sage: padic_printing.mode('series')
        """
        if sep is None:
            return self._sep
        else:
            self._sep = str(sep)

    def alphabet(self, alphabet = None):
        """
        EXAMPLES:
            sage: padic_printing.alphabet("abc")
            sage: padic_printing.mode('digits')
            sage: repr(Qp(3)(1234))
            '...bcaacab'

            sage: padic_printing.mode('series')
        """
        if alphabet is None:
            return self._alphabet
        else:
            self._alphabet = list(alphabet)

_printer_defaults = pAdicPrinterDefaults()


cdef class pAdicPrinter_class(SageObject):
    def __init__(self, ring, mode, pos, pname, unram_name, var_name, max_ram_terms, max_unram_terms, max_terse_terms, sep, alphabet):
        global _printer_defaults
        self.ring = ring
        self.prime_pow = ring.prime_pow
        from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
        self.base = isinstance(ring, pAdicBaseGeneric)
        if alphabet is None:
            self.alphabet = _printer_defaults._alphabet
        else:
            self.alphabet = alphabet
        # note that self.pos is reset to True if mode == 'digits'
        if pos is None:
            self.pos = _printer_defaults._pos
        else:
            self.pos = pos
        if mode is None:
            mode = _printer_defaults._mode
        if mode == 'val-unit':
            self.mode = val_unit
        elif mode == 'series':
            self.mode = series
        elif mode == 'terse':
            self.mode = terse
        elif mode == 'digits':
            if len(self.alphabet) < self.prime_pow.prime or (not self.base and ring.inertia_degree() != 1):
                self.mode = bars
            else:
                self.mode = digits
                self.pos = True
        elif mode == 'bars':
            self.mode = bars
        else:
            raise ValueError, "printing mode must be one of 'val-unit', 'series', 'terse', 'digits' or 'bars'"
        if pname is None:
            self.pname = ring._uniformizer_print()
        else:
            self.pname = pname
        if unram_name is None:
            self.unram_name = ring._unram_print()
        else:
            self.unram_name = unram_name
        if var_name is None:
            self.var_name = ring.variable_name()
        else:
            self.var_name = var_name
        if sep is None:
            self.sep = _printer_defaults._sep
        else:
            self.sep = sep
        if max_ram_terms is not None:
            self.max_ram_terms = max_ram_terms
            if self.max_ram_terms < 0:
                raise ValueError, "max_ram_terms must be positive and fit in a long"
        else:
            self.max_ram_terms = _printer_defaults._max_ram_terms
        if max_unram_terms is not None:
            self.max_unram_terms = max_unram_terms
            if self.max_unram_terms < 0:
                raise ValueError, "max_unram_terms must be positive and fit in a long"
        else:
            self.max_unram_terms = _printer_defaults._max_unram_terms
        if max_terse_terms is not None:
            self.max_terse_terms = max_terse_terms
            if self.max_terse_terms < 0:
                raise ValueError, "max_terse_terms must be positive and fit in a long"
        else:
            self.max_terse_terms = _printer_defaults._max_terse_terms

    def __reduce__(self):
        return pAdicPrinter, (self.ring, \
                              self._print_mode(), \
                              self.pos, \
                              self.pname, \
                              self.unram_name, \
                              self.var_name, \
                              None if self.max_ram_terms == -1 else self.max_ram_terms, \
                              None if self.max_unram_terms == -1 else self.max_unram_terms, \
                              None if self.max_terse_terms == -1 else self.max_terse_terms, \
                              self.sep, \
                              self.alphabet)

    def _repr_(self):
        return "%s printer for %s"%(self._print_mode(), self.ring)

    def __enter__(self):
        self.old = self.ring._printer
        self.ring._printer = self

    def __exit__(self, type, value, traceback):
        self.ring._printer = self.old

    def _pos(self):
        return self.pos

    def _set_pos(self, pos):
        self.pos = pos

    def _ring(self):
        return self.ring

    def _uniformizer_name(self):
        return self.pname

    def _print_mode(self):
        if self.mode == val_unit:
            return 'val-unit'
        elif self.mode == series:
            return 'series'
        elif self.mode == terse:
            return 'terse'
        elif self.mode == digits:
            return 'digits'
        elif self.mode == bars:
            return 'bars'

    def _base_p_list(self, value, pos):
        cdef Integer _value = Integer(value)
        return self.base_p_list(_value.value, pos)

    cdef base_p_list(self, mpz_t value, bint pos):
        cdef mpz_t tmp, halfp
        cdef int neg, curpower
        cdef Integer list_elt
        cdef unsigned long preccap = self.prime_pow.prec_cap
        ans = PyList_New(0)
        mpz_init_set(tmp, value)


        list_elt = PY_NEW(Integer)
        mpz_set(list_elt.value, value)
        if pos:
            while mpz_sgn(tmp) != 0:
                list_elt = PY_NEW(Integer)
                mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
                mpz_sub(tmp, tmp, list_elt.value)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                PyList_Append(ans, list_elt)
        else:
            neg = 0
            curpower = preccap
            mpz_init(halfp)
            mpz_fdiv_q_2exp(halfp, self.prime_pow.prime.value, 1)
            while mpz_sgn(tmp) != 0:
                curpower -= 1
                list_elt = PY_NEW(Integer)
                mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
                if mpz_cmp(list_elt.value, halfp) > 0:
                    mpz_sub(list_elt.value, list_elt.value, self.prime_pow.prime.value)
                    neg = 1
                else:
                    neg = 0
                mpz_sub(tmp, tmp, list_elt.value)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                if neg == 1:
                    if mpz_cmp(tmp, self.prime_pow.pow_mpz_t_tmp(curpower)[0]) >= 0:
                        mpz_sub(tmp, tmp, self.prime_pow.pow_mpz_t_tmp(curpower)[0])
                PyList_Append(ans, list_elt)
            mpz_clear(halfp)
        mpz_clear(tmp)
        return ans

    def repr_gen(self, elt, do_latex, pos = None, mode = None, pname = None):
        cdef int _mode
        cdef bint _pos
        if mode is None:
            _mode = self.mode
        elif mode == 'val-unit':
            _mode = val_unit
        elif mode == 'series':
            _mode = series
        elif mode == 'terse':
            _mode = terse
        elif mode == 'digits':
            _mode = digits
        elif mode == 'bars':
            _mode = bars
        else:
            raise ValueError, "printing mode must be one of 'val-unit', 'series', 'terse', 'bars', or 'digits'"
        if pos is None:
            _pos = self.pos
        else:
            _pos = pos
        if pname is None:
            pprint = self.pname
        else:
            pprint = str(pname)
        return self._repr_gen(elt, do_latex, _pos, _mode, pprint)

    cdef _repr_gen(self, pAdicGenericElement elt, bint do_latex, bint pos, int mode, pname):
        r"""
        Prints a string representation of the element.  See __init__ for more details on print modes.

        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel','val-unit'); a = R(364); a
            7 * 52 + O(7^5)
            sage: print a.str('terse')
            364 + O(7^5)
            sage: print a.str('series')
            3*7 + 7^3 + O(7^5)
            sage: K = Qp(7,4,'capped-rel','val-unit'); a = K(364); a
            7 * 52 + O(7^5)
            sage: print a.str('series')
            3*7 + 7^3 + O(7^5)
        """
        cdef Py_ssize_t i
        if elt._is_exact_zero():
            return "0"
        if elt._is_inexact_zero():
            if mode == val_unit or mode == series:
                s = "O(%s"%(pname)
            elif mode == terse:
                s = "0 + O(%s"%(pname)
            else: # mode == digits or bars
                s = "..."
        elif mode == val_unit:
            if do_latex:
                if elt.valuation() == 0:
                    s = "%s + O(%s"%(self._repr_spec(elt, do_latex, pos, terse, 0, pname), pname)
                elif elt.valuation() == 1:
                    s = "%s \\cdot %s + O(%s"%(pname, self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, pname), pname)
                else:
                    s = "%s^{%s} \\cdot %s + O(%s"%(pname, elt.valuation(), self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, pname), pname)
            else:
                if elt.valuation() == 0:
                    s = "%s + O(%s"%(self._repr_spec(elt, do_latex, pos, terse, 0, pname), pname)
                elif elt.valuation() == 1:
                    s = "%s * %s + O(%s"%(pname, self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, pname), pname)
                else:
                    s = "%s^%s * %s + O(%s"%(pname, elt.valuation(), self._repr_spec(elt.unit_part(), do_latex, pos, terse, 1, pname), pname)
        elif mode == digits:
            if self.base:
                L = self.base_p_list((<Integer>elt.unit_part().lift()).value, True)
            else:
                L = elt._ext_p_list(self.pos)
            if self.max_ram_terms != -1:
                L = L[:self.max_ram_terms]
            L.reverse()
            # The following step should work since mode is only allowed to be digits in the case of totally ramified extensions
            # with primes smaller than the length of the alphabet
            L = [self.alphabet[a] for a in L]
            n = elt.valuation()
            if n > 0:
                L += [self.alphabet[0]]*n
            elif n < 0:
                print L
                print n
                print L[:n]
                print L[n:]
                L = L[:n] + ['.'] + L[n:]
            s = "".join(L)
            s = "..." + s
        elif mode == bars:
            if self.base:
                L = self.base_p_list((<Integer>elt.unit_part().lift()).value, self.pos)
            else:
                L = elt._ext_p_list(self.pos)
            if self.max_ram_terms != -1:
                L = L[:self.max_ram_terms]
            L.reverse()
            L = [str(a) for a in L]
            n = elt.valuation()
            if n > 0:
                L += ['0']*n
            elif n < 0:
                L = L[:n] + ['.'] + L[n:]
            s = self.sep.join(L)
            s = "..." + s
        else: # mode == terse or series
            s = "%s + O(%s"%(self._repr_spec(elt, do_latex, pos, mode, 0, pname), pname)
        if mode != bars and mode != digits:
            if elt.precision_absolute() == 1:
                s += ")"
            else:
                if do_latex:
                    s += "^{%s})"%(elt.precision_absolute())
                else:
                    s += "^%s)"%(elt.precision_absolute())
        return s

    cdef _repr_spec(self, pAdicGenericElement elt, bint do_latex, bint pos, int mode, bint paren, pname):
        """
        Should not be called if elt is an exact or inexact zero
        """
        cdef Integer lift_z, pprec
        cdef int ZZ_pEX
        cdef Py_ssize_t i, j
        cdef long val
        #cdef bint ellipsis = 0
        cdef ellipsis_unram
        if self.base:
            if mode == terse:
                if elt.valuation() >= 0:
                    lift_z = <Integer> elt.lift()
                    if not pos:
                        pprec = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>elt.precision_absolute()).value))
                        if lift_z > pprec / 2:
                            mpz_sub(lift_z.value, lift_z.value, pprec.value)
                            if paren:
                                return "(%s)"%(lift_z)
                    return "%s"%(lift_z)
                else:
                    ppow1 = str(self.prime_pow.pow_Integer(-mpz_get_si((<Integer>elt.valuation()).value)))
                    if do_latex:
                        ppow2 = "%s^{%s}"%(pname, -elt.valuation())
                    else:
                        ppow2 = "%s^%s"%(pname, -elt.valuation())
                    if len(ppow1) < len(ppow2):
                        ppow = ppow1
                    else:
                        ppow = ppow2
                    if do_latex:
                        return "\frac{%s}{%s}"%(self._repr_spec(elt.unit_part(), do_latex, pos, terse, 0, pname), ppow)
                    elif paren:
                        return "(%s/%s)"%(self._repr_spec(elt.unit_part(), do_latex, pos, terse, 0, pname), ppow)
                    else:
                        return "%s/%s"%(self._repr_spec(elt.unit_part(), do_latex, pos, terse, 0, pname), ppow)
            else: # mode == series
                slist = self.base_p_list((<Integer>elt.unit_part().lift()).value, pos)
                slist, ellipsis = self._truncate_list(slist, self.max_ram_terms, 0)
                s = ""
                exp = elt.valuation()
                for a in slist:
                    if a != 0:
                        if a < 0:
                            if len(s) == 0:
                                s = "-"
                            else:
                                s += " - "
                            a = -a
                        elif len(s) != 0:
                            s += " + "
                        s += self._co_dot_var(a, pname, exp, do_latex)
                    exp += 1
                if ellipsis:
                    s += self._plus_ellipsis(do_latex)
        else: # not self.base
            if mode == terse:
                if elt.parent().is_capped_relative():
                    poly, k = elt._ntl_rep_abs()
                    s = repr(poly)
                else:
                    s = repr(elt._ntl_rep())
                L = s.split("] [") # this splits a ZZ_pEX into the ZZ_pX components
                ZZ_pEX = L[0].count("[") # will equal 2 if elt was a ZZ_pEX element, 1 if it was a ZZ_pX element
                L[0] = L[0].replace("[","")
                L[-1] = L[-1].replace("]","")
                if ZZ_pEX == 2:
                    L = [a.split() for a in L]
                    L = [[("" if b == "0" else b) for b in a] for a in L]
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, "")
                    raise NotImplementedError
                else:
                    L = L[0].split()
                    L = [("" if b == "0" else b) for b in L]
                    L, ellipsis = self._truncate_list(L, self.max_terse_terms, "")
                    s = ""
                    pn = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>elt.precision_absolute()).value))
                    if elt.parent().is_capped_relative():
                        pk = self.prime_pow.pow_Integer(mpz_get_ui((<Integer>k).value))
                    else:
                        pk = Integer(1)
                    for i from 0 <= i < len(L):
                        if L[i] != "":
                            a = Integer(L[i])
                            if not pos and 2*a > pn:
                                a = (pn - a) / pk
                                if s == "":
                                    s = "-%s"%(a)
                                elif a == 1:
                                    s += " - "
                                else:
                                    s += " - %s"%(a)
                                s += self._dot_var(self.var_name, i, do_latex)
                            elif a == pk:
                                if s != "":
                                    s += " + "
                                s += self._var(self.var_name, i, do_latex)
                            else:
                                a = a / pk
                                if s == "":
                                    s = "%s"%a
                                else:
                                    s += " + %s"%a
                                s += self._dot_var(self.var_name, i, do_latex)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
            else: # series
                s = ""
                L = elt._ext_p_list(pos)
                val = elt.valuation_c()
                # since elt was not supposed to be zero, this should give a non-empty list.
                if len(L) == 0:
                    raise RuntimeError, "repr_spec called on zero"
                if isinstance(L[0], list): # unramified part to the extension
                    if self.unram_name is None:
                        raise RuntimeError, "need to have specified a name for the unramified variable"
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, [])
                    for i from 0 <= i < len(L):
                        L[i], ellipsis_unram = self._truncate_list(L[i], self.max_unram_terms, 0)
                        term = self._print_list_as_poly(L[i], do_latex, self.unram_name, 0, 0)
                        if len(term) > 0:
                            if ellipsis_unram:
                                term += self._plus_ellipsis(do_latex)
                            exp = i + val
                            if (not do_latex and term.find(" ") != -1) or (do_latex and (term.find(" + ") != -1 or term.find(" - ") != -1)):
                                if len(s) > 0:
                                    s += " + (" + term + ")"
                                else:
                                    s = "(" + term + ")"
                                s += self._dot_var(pname, exp, do_latex)
                            else: # since len(term) > 0, len(L[i]) == 1.
                                if term == "1":
                                    if len(s) > 0:
                                        s += " + "
                                    s += self._var(pname, exp, do_latex)
                                    continue
                                elif term == "-1":
                                    if len(s) > 0:
                                        s += " - "
                                    else:
                                        s = "-"
                                    s += self._var(pname, exp, do_latex)
                                    continue
                                elif len(s) > 0:
                                    if term[0] == "-":
                                        s += " - " + term[1:]
                                    else:
                                        s += " + " + term
                                    s += self._dot_var(pname, exp, do_latex)
                                else:
                                    s = term
                                    s += self._dot_var(pname, exp, do_latex)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
                else: # L[0] is not a list, so no unramified printing required
                    L, ellipsis = self._truncate_list(L, self.max_ram_terms, 0)
                    s = self._print_list_as_poly(L, do_latex, pname, val, 1)
                    if ellipsis:
                        s += self._plus_ellipsis(do_latex)
        if paren and s.find(" ") != -1:
            s = "(" + s + ")"
        return s

    cdef _var(self, x, exp, do_latex):
        if exp == 0:
            return "1"
        if exp == 1:
            return str(x)
        if do_latex:
            return "%s^{%s}"%(x, exp)
        else:
            return "%s^%s"%(x, exp)

    cdef _dot_var(self, x, exp, do_latex):
        if exp == 0:
            return ""
        if exp == 1:
            if do_latex:
                return " \\cdot %s"%(x)
            else:
                return "*%s"%(x)
        if do_latex:
            return " \\cdot %s^{%s}"%(x, exp)
        else:
            return "*%s^%s"%(x, exp)

    cdef _co_dot_var(self, co, x, exp, do_latex):
        """
        co should be greater than 0
        """
        if exp == 0:
            return "%s"%co
        if exp == 1:
            if co == 1:
                return "%s"%x
            if do_latex:
                return "%s \\cdot %s"%(co, x)
            else:
                return "%s*%s"%(co, x)
        if co == 1:
            if do_latex:
                return "%s^{%s}"%(x, exp)
            else:
                return "%s^%s"%(x, exp)
        if do_latex:
            return "%s \\cdot %s^{%s}"%(co, x, exp)
        else:
            return "%s*%s^%s"%(co, x, exp)

    cdef _plus_ellipsis(self, do_latex):
        if do_latex:
            return " + \\cdots"
        else:
            return " + ..."

    cdef _truncate_list(self, L, max_terms, zero):
        cdef bint ellipsis = 0
        if max_terms == -1:
            return L, ellipsis
        if len(L) == 0:
            return L, ellipsis
        cdef Py_ssize_t i, nonzero_index
        cdef Py_ssize_t count = 0
        for i from 0 <= i < len(L):
            if L[i] != zero:
                count += 1
                if count > max_terms:
                    ellipsis = 1
                    return L[:nonzero_index+1], ellipsis
                nonzero_index = i
        return L, ellipsis

    cdef _print_list_as_poly(self, L, bint do_latex, polyname, long expshift, bint increasing):
        s = ""
        cdef Py_ssize_t j
        cdef long exp
        if increasing:
            for j from 0 <= j < len(L):
                exp = j + expshift
                s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        else:
            for j from len(L) - 1 >= j >= 0:
                exp = j + expshift
                s = self._print_term_of_poly(s, L[j], do_latex, polyname, exp)
        return s

    cdef _print_term_of_poly(self, s, coeff, bint do_latex, polyname, long exp):
        if coeff < 0:
            if len(s) > 0:
                s += " - "
            else:
                s = "-"
            coeff = -coeff
            s += self._co_dot_var(coeff, polyname, exp, do_latex)
        elif coeff > 0:
            if len(s) > 0:
                s += " + "
            s += self._co_dot_var(coeff, polyname, exp, do_latex)
        return s
