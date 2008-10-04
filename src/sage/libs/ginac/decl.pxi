###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

# NOTE: Because of the except+'s below, i.e., C++ exception handling,
# we do *not* have to use _sig_on and _sig_off. We do use it a little
# in the actual pyx code to catch control-c for long running functions.

cdef extern from "ginac_wrap.h":
    void ginac_pyinit_Integer(object)
    void ginac_pyinit_Float(object)
    void ginac_pyinit_I(object)

    ctypedef struct GBasic "basic":
        unsigned int gethash()
        int compare(GBasic other)

    ctypedef struct GSymbol "symbol":
        pass

    GSymbol* GSymbol_construct_str "Construct_p<symbol, char*>" \
            (void *mem, char* m)

    void GSymbol_destruct "Destruct<symbol>"(GSymbol *mem)

    object GSymbol_to_str "_to_PyString<symbol>"(GSymbol *s)

    ctypedef struct GEx "ex":
        unsigned int gethash()        except +
        int compare(GEx other)        except +
        GEx expand(unsigned int opt)  except +
        GEx collect(GEx s, bint dist) except +
        GEx diff(GSymbol s, int d)    except +
        GEx series(GEx s, int order, unsigned options) except +
        bint is_zero()                except +
        bint match(GEx pattern)       except +
        bint has(GEx pattern)         except +
        GEx subs(GEx expr)            except +
        GEx coeff(GEx expr, int n)    except +
        GEx lcoeff(GEx expr)          except +
        GEx tcoeff(GEx expr)          except +
        int degree(GEx expr)          except +
        int ldegree(GEx expr)         except +
        GEx rhs()                     except +
        GEx lhs()                     except +
        int nops()                    except +
        GEx op(int i)                 except +
        GEx eval(int level)           except +
        GEx evalf(int level)          except +
        GEx conjugate()               except +
        GEx real_part()               except +
        GEx imag_part()               except +

    # Numericals
    bint is_a_numeric "is_a<numeric>" (GEx e)
    # given a GEx that is known to be a numeric, return reference to
    # the underlying PyObject*.
    object py_object_from_numeric(GEx e)     except +

    # Algorithms
    GEx g_gcd "gcd"(GEx a, GEx b) except +

    # Pattern matching wildcards
    GEx g_wild "wild"(unsigned int label) except +

    # Series back to poly
    GEx series_to_poly(GEx e) except +
    bint is_a_series "is_a<pseries>" (GEx e)

    # Relations
    ctypedef enum operators "relational::operators":
        equal, not_equal, less, less_or_equal, greater, greater_or_equal
    bint is_negative(GEx x)                  except +
    bint is_a_relational "is_a<relational>" (GEx e)
    bint relational_to_bool(GEx e)
    operators relational_operator(GEx e)
    operators switch_operator(operators op)
    GEx relational(GEx lhs, GEx rhs, operators o)
    GEx g_lt "LT_WRAP" (GEx left, GEx right) except +
    GEx g_eq "EQ_WRAP" (GEx left, GEx right) except +
    GEx g_gt "GT_WRAP" (GEx left, GEx right) except +
    GEx g_le "LE_WRAP" (GEx left, GEx right) except +
    GEx g_ne "NE_WRAP" (GEx left, GEx right) except +
    GEx g_ge "GE_WRAP" (GEx left, GEx right) except +

    # Constants
    GEx g_Pi "Pi"
    GEx g_Catalan "Catalan"
    GEx g_Euler "Euler"

    # Destructor and constructor
    void GEx_destruct "Destruct<ex>"(GEx *mem) except +
    GEx* GEx_construct_symbol "Construct_p<ex, symbol>" \
            (void *mem, GSymbol m) except +
    GEx* GEx_construct_ex "Construct_p<ex, ex>" (void *mem, GEx m) except +
    GEx* GEx_construct_long "Construct_p<ex, long>" (void *mem, long n) except +
    GEx* GEx_construct_double "Construct_p<ex, double>" \
            (void *mem, double d) except +

    GEx* GEx_construct_pyobject "ASSIGN_WRAP" (GEx mem, object n)

    # Conversions
    double GEx_to_double(GEx e, int* success) except +
    object GEx_to_str "_to_PyString<ex>"(GEx *s)

    bint is_a_symbol "is_a<symbol>" (GEx e)
    GSymbol ex_to_symbol "ex_to<symbol>" (GEx e)

    ctypedef struct GExVector "exvector":
        void push_back(GEx)
        int size()
        GEx at(int i)

    ctypedef struct GExSetIter "std::set<ex, ex_is_less>::const_iterator":
        void inc "operator++" ()
        GEx obj "operator*" ()
        bint is_not_equal "operator!=" (GExSetIter i)

    ctypedef struct GExSet "std::set<ex, ex_is_less>":
        GExSetIter begin()
        GExSetIter end()

    void g_list_symbols "list_symbols" (GEx e, GExSet s)

    # Arithmetic
    int ginac_error()
    GEx gadd "ADD_WRAP" (GEx left, GEx right) except +
    GEx gsub "SUB_WRAP" (GEx left, GEx right) except +
    GEx gmul "MUL_WRAP" (GEx left, GEx right) except +
    GEx gdiv "DIV_WRAP" (GEx left, GEx right) except +
    GEx g_pow "pow" (GEx left, GEx exp)      except +

    GSymbol get_symbol(char* s)              except +

    GEx g_abs "GiNaC::abs" (GEx x)           except +
    GEx g_step "GiNaC::step" (GEx x)	     except +  # step function
    GEx g_csgn "GiNaC::csgn" (GEx x)	     except + # complex sign
    GEx g_conjugate "GiNaC::conjugate" (GEx x)	except + # complex conjugation
    GEx g_real_part "GiNaC::real_part" (GEx x)	except + # real part
    GEx g_imag_part "GiNaC::imag_part" (GEx x)	except + # imaginary part
    GEx g_sqrt "GiNaC::sqrt" (GEx x)	except +  # square root (not a GiNaC function, rather an alias for pow(x, numeric(1, 2)))
    GEx g_sin "GiNaC::sin" (GEx x)	except + # sine
    GEx g_cos "GiNaC::cos" (GEx x)	except + # cosine
    GEx g_tan "GiNaC::tan" (GEx x)	except + # tangent
    GEx g_asin "GiNaC::asin" (GEx x)	except + # inverse sine
    GEx g_acos "GiNaC::acos" (GEx x)	except + # inverse cosine
    GEx g_atan "GiNaC::atan" (GEx x)	except + # inverse tangent
    GEx g_atan2 "GiNaC::atan2" (GEx y, GEx x) except + 	# inverse tangent with two arguments
    GEx g_sinh "GiNaC::sinh" (GEx x)	except + # hyperbolic sine
    GEx g_cosh "GiNaC::cosh" (GEx x)	except + # hyperbolic cosine
    GEx g_tanh "GiNaC::tanh" (GEx x)	except + # hyperbolic tangent
    GEx g_asinh "GiNaC::asinh" (GEx x)	except + # inverse hyperbolic sine
    GEx g_acosh "GiNaC::acosh" (GEx x)	except + # inverse hyperbolic cosine
    GEx g_atanh "GiNaC::atanh" (GEx x)	except + # inverse hyperbolic tangent
    GEx g_exp "GiNaC::exp" (GEx x)	except + # exponential function
    GEx g_log "GiNaC::log" (GEx x)	except + # natural logarithm
    GEx g_Li2 "GiNaC::Li2" (GEx x)	except + # dilogarithm
    GEx g_Li "GiNaC::Li" (GEx m, GEx x)	except + # classical polylogarithm as well as multiple polylogarithm
    GEx g_G "GiNaC::G" (GEx a, GEx y)	except + # multiple polylogarithm
    GEx g_G2 "GiNaC::G" (GEx a, GEx s, GEx y)	except + # multiple polylogarithm with explicit signs for the imaginary parts
    GEx g_S "GiNaC::S" (GEx n, GEx p, GEx x)	except + # Nielsen's generalized polylogarithm
    GEx g_H "GiNaC::H" (GEx m, GEx x)	        except + # harmonic polylogarithm
    GEx g_zeta "GiNaC::zeta" (GEx m)	        except + # Riemann's zeta function as well as multiple zeta value
    GEx g_zeta2 "GiNaC::zeta" (GEx m, GEx s)	except + # alternating Euler sum
    GEx g_zetaderiv "GiNaC::zetaderiv" (GEx n, GEx x)	except + # derivatives of Riemann's zeta function
    GEx g_tgamma "GiNaC::tgamma" (GEx x)	except + # gamma function
    GEx g_lgamma "GiNaC::lgamma" (GEx x)	except + # logarithm of gamma function
    GEx g_beta "GiNaC::beta" (GEx x, GEx y)	except + # beta function (tgamma*tgamma(y)/tgamma(x+y))
    GEx g_psi "GiNaC::psi" (GEx x)	        except + # psi (digamma) function
    GEx g_psi2 "GiNaC::psi" (GEx n, GEx x)	except + # derivatives of psi function (polygamma functions)
    GEx g_factorial "GiNaC::factorial" (GEx n)	except + # factorial function n!
    GEx g_binomial "GiNaC::binomial" (GEx n, GEx k)	except + # binomial coefficients
    GEx g_Order "GiNaC::Order" (GEx x)	        except + # order term function in truncated power series


    ctypedef struct GFunction "function":
        pass

    GEx g_function_evalv(unsigned int serial, GExVector) except +
    GEx g_function_eval0(unsigned int serial) except +
    GEx g_function_eval1(unsigned int serial, GEx) except +
    GEx g_function_eval2(unsigned int serial, GEx, GEx) except +
    GEx g_function_eval3(unsigned int serial, GEx, GEx, GEx) except +
    GEx g_function_eval4(unsigned int serial, GEx, GEx, GEx, GEx) except +
    GEx g_function_eval5(unsigned int serial, GEx, GEx, GEx, GEx, GEx) except +
    GEx g_function_eval6(unsigned int serial, GEx, GEx, GEx, GEx, GEx, \
            GEx) except +
    GEx g_function_eval7(unsigned int serial, GEx, GEx, GEx, GEx, GEx, \
            GEx, GEx) except +
    GEx g_function_eval8(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx) except +
    GEx g_function_eval9(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx) except +
    GEx g_function_eval10(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx, GEx) except +
    GEx g_function_eval11(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx, GEx, GEx) except +
    GEx g_function_eval12(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx, GEx, GEx, GEx) except +
    GEx g_function_eval13(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx, GEx, GEx, GEx, GEx) except +
    GEx g_function_eval14(unsigned int serial, GEx, GEx, GEx, GEx, GEx, GEx,
            GEx, GEx, GEx, GEx, GEx, GEx, GEx, GEx) except +

    ctypedef struct GFunctionOpt "function_options":
        void set_python_func()
        GFunctionOpt eval_func(object f)
        GFunctionOpt evalf_func(object f)
        GFunctionOpt conjugate_func(object f)
        GFunctionOpt real_part_func(object f)
        GFunctionOpt imag_part_func(object f)
        GFunctionOpt derivative_func(object f)
        GFunctionOpt power_func(object f)
        GFunctionOpt series_func(object f)
        GFunctionOpt print_func(object f)

    ctypedef struct GFunctionOptVector "vector<function_options>":
        int size()
        GFunctionOpt index "operator[]" (int ind)

    void g_foptions_assign "ASSIGN_WRAP" (GFunctionOpt, GFunctionOpt)

    GFunctionOpt g_function_options "GiNaC::function_options" \
            (char *m)
    GFunctionOpt g_function_options_args "GiNaC::function_options" \
            (char *m, unsigned nargs)
    unsigned g_register_new "GiNaC::function::register_new" (GFunctionOpt opt)

    unsigned find_function "GiNaC::function::find_function" (char* name,
            unsigned nargs) except +ValueError

    GFunctionOptVector g_registered_functions \
            "GiNaC::function::registered_functions" ()

    unsigned cos_serial "cos_SERIAL::serial"
