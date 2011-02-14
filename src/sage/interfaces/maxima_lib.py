r"""
Interface to Maxima

Maxima is a free GPL'd general purpose computer algebra system
whose development started in 1968 at MIT. It contains symbolic
manipulation algorithms, as well as implementations of special
functions, including elliptic functions and generalized
hypergeometric functions. Moreover, Maxima has implementations of
many functions relating to the invariant theory of the symmetric
group `S_n`. (However, the commands for group invariants,
and the corresponding Maxima documentation, are in French.) For many
links to Maxima documentation see
http://maxima.sourceforge.net/docs.shtml/.

AUTHORS:

- William Stein (2005-12): Initial version

- David Joyner: Improved documentation

- William Stein (2006-01-08): Fixed bug in parsing

- William Stein (2006-02-22): comparisons (following suggestion of
  David Joyner)

- William Stein (2006-02-24): *greatly* improved robustness by adding
  sequence numbers to IO bracketing in _eval_line

If the string "error" (case insensitive) occurs in the output of
anything from Maxima, a RuntimeError exception is raised.
"""

#*****************************************************************************
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

from sage.symbolic.ring import SR

from sage.libs.ecl import *

from maxima_abstract import MaximaAbstract, MaximaAbstractFunction, MaximaAbstractElement, MaximaAbstractFunctionElement, MaximaAbstractElementFunction

## We begin here by initializing maxima in library mode
ecl_eval("(setf *load-verbose* NIL)")
ecl_eval("(require 'maxima)")
ecl_eval("(in-package :maxima)")
ecl_eval("(setq $nolabels t))")
ecl_eval("(defvar *MAXIMA-LANG-SUBDIR* NIL)")
ecl_eval("(set-locale-subdir)")
ecl_eval("(set-pathnames)")
ecl_eval("(defun add-lineinfo (x) x)")
#the following is a direct adaption of the definition of "retrieve" in the Maxima file
#macsys.lisp. This routine is normally responsible for displaying a question and
#returning the answer. We change it to throw an error in which the text of the question
#is included. We do this by running exactly the same code as in the original definition
#of "retrieve", but with *standard-output* redirected to a string.
ecl_eval(r"""
(defun retrieve (msg flag &aux (print? nil))
  (declare (special msg flag print?))
  (or (eq flag 'noprint) (setq print? t))
  (error
      (concatenate 'string "Maxima asks: "
      (string-trim '(#\Newline)
      (with-output-to-string (*standard-output*)
      (cond ((not print?)
             (setq print? t)
             (princ *prompt-prefix*)
             (princ *prompt-suffix*)
             )
            ((null msg)
             (princ *prompt-prefix*)
             (princ *prompt-suffix*)
             )
            ((atom msg)
             (format t "~a~a~a" *prompt-prefix* msg *prompt-suffix*)
             )
            ((eq flag t)
             (princ *prompt-prefix*)
             (mapc #'princ (cdr msg))
             (princ *prompt-suffix*)
             )
            (t
             (princ *prompt-prefix*)
             (displa msg)
             (princ *prompt-suffix*)
             )
      ))))
  )
)
""")

ecl_eval('(defparameter *dev-null* (make-two-way-stream (make-concatenated-stream) (make-broadcast-stream)))')
ecl_eval('(defun principal nil (error "Divergent Integral"))')
ecl_eval("(setf $errormsg nil)")

#ecl_eval(r"(defun tex-derivative (x l r) (tex (if $derivabbrev (tex-dabbrev x) (tex-d x '\partial)) l r lop rop ))")

#ecl_eval('(defun ask-evod (x even-odd)(error "Maxima asks a question"))')
#ecl_eval('(defun ask-integerp (x)(error "Maxima asks a question"))')
#ecl_eval('(defun ask-declare (x property)(error "Maxima asks a question"))')
#ecl_eval('(defun ask-prop (object property fun-or-number)(error "Maxima asks a question"))')
#ecl_eval('(defun asksign01 (a)(error "Maxima asks a question"))')
#ecl_eval('(defun asksign (x)(error "Maxima asks a question"))')
#ecl_eval('(defun asksign1 ($askexp)(error "Maxima asks a question"))')
#ecl_eval('(defun ask-greateq (x y)(error "Maxima asks a question"))')
#ecl_eval('(defun askinver (a)(error "Maxima asks a question"))')
#ecl_eval('(defun npask (exp)(error "Maxima asks a question"))')

ecl_eval("(setf original-standard-output *standard-output*)")
ecl_eval("(setf *standard-output* *dev-null*)")
#ecl_eval("(setf *error-output* *dev-null*)")

# display2d -- no ascii art output
# keepfloat -- don't automatically convert floats to rationals
init_code = ['display2d : false', 'domain : complex', 'keepfloat : true', 'load(to_poly_solver)', 'load(simplify_sum)']
# Turn off the prompt labels, since computing them *very
# dramatically* slows down the maxima interpret after a while.
# See the function makelabel in suprv1.lisp.
# Many thanks to andrej.vodopivec@gmail.com and also
# Robert Dodier for figuring this out!
# See trac # 6818.
init_code.append('nolabels : true')
for l in init_code:
    ecl_eval("#$%s$"%l)
## To get more debug information uncomment the next line
## should allow to do this through a method
#ecl_eval("(setf *standard-output* original-standard-output)")

# This returns an EclObject
maxima_eval=ecl_eval("""
(defun maxima-eval( form )
    (let ((result (catch 'macsyma-quit (cons 'maxima_eval (meval form)))))
        ;(princ (list "result=" result))
        ;(terpri)
        ;(princ (list "$error=" $error))
        ;(terpri)
        (cond
            ((and (consp result) (eq (car result) 'maxima_eval)) (cdr result))
            ((eq result 'maxima-error)
                (let ((the-jig (process-error-argl (cddr $error))))
                    (mapc #'set (car the-jig) (cadr the-jig))
                    (error (concatenate 'string "Error executing code in Maxima: "
                       (with-output-to-string (stream)
                           (apply #'mformat stream (cadr $error) (caddr the-jig)))))
                ))
            (t
                (let ((the-jig (process-error-argl (cddr $error))))
                    (mapc #'set (car the-jig) (cadr the-jig))
                    (error (concatenate 'string "Maxima condition. result:" (princ-to-string result) "$error:"
                       (with-output-to-string (stream)
                           (apply #'mformat stream (cadr $error) (caddr the-jig)))))
                ))
        )
    )
)
""")

maxima_lib_instances = 0

# The Maxima string function can change the structure of its input
#maxprint=EclObject("$STRING")
maxprint=EclObject("(defun mstring-for-sage (form) (coerce (mstring form) 'string))").eval()
meval=EclObject("MEVAL")
msetq=EclObject("MSETQ")
mlist=EclObject("MLIST")
mequal=EclObject("MEQUAL")
cadadr=EclObject("CADADR")

max_integrate=EclObject("$INTEGRATE")
max_sum=EclObject("$SUM")
max_simplify_sum=EclObject("$SIMPLIFY_SUM")
max_ratsimp=EclObject("$RATSIMP")
max_limit=EclObject("$LIMIT")
max_tlimit=EclObject("$TLIMIT")
max_plus=EclObject("$PLUS")
max_minus=EclObject("$MINUS")
max_use_grobner=EclObject("$USE_GROBNER")
max_to_poly_solve=EclObject("$TO_POLY_SOLVE")

def stdout_to_string(s):
    return ecl_eval("(with-output-to-string (*standard-output*) (maxima-eval #$%s$))"%s).python()[1:-1]

def max_to_string(s):
     return maxprint(s).python()[1:-1]

my_mread=ecl_eval("""
(defun my-mread (cmd)
  (caddr (mread (make-string-input-stream cmd))))
""")

def parse_max_string(l):
  return my_mread('"%s;"'%l)

class MaximaLib(MaximaAbstract):
    """
    Interface to Maxima as a Library.
    """
    def __init__(self):
        """
        Create an instance of the Maxima interpreter.

        TESTS::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib == loads(dumps(maxima_lib))
            True

        We make sure labels are turned off (see trac 6816)::

            sage: 'nolabels : true' in maxima_lib._MaximaLib__init_code
            True
        """
        global maxima_lib_instances
        if maxima_lib_instances > 0:
            raise RuntimeError, "Maxima interface in library mode can only be instantiated once"
        maxima_lib_instances += 1

        global init_code
        self.__init_code = init_code

        ## The name should definitely be changed to maxima_lib, however much more changes are then needed elsewhere
        ## With maxima, more things are fine, but for example _maxima_init_ gets called in calculus.calculus and the classic interface gets initialized (not started, it is already initialized by default, so that is not really a big deal)
        MaximaAbstract.__init__(self,"maxima_lib")
        self.__seq = 0

    def _coerce_from_special_method(self, x):
        if isinstance(x, EclObject):
            return MaximaLibElement(self,self._create(x))
        else:
            return MaximaAbstract._coerce_from_special_method(self,x)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.__reduce__()
            (<function reduce_load_MaximaLib at 0x...>, ())
        """
        return reduce_load_MaximaLib, tuple([])

    # This outputs a string
    def eval(self, line, locals=None, reformat=True, **kwds):
        result = ''
        while line:
            ind_dollar=line.find("$")
            ind_semi=line.find(";")
            if ind_dollar == -1 or (ind_semi >=0 and ind_dollar > ind_semi):
                if ind_semi == -1:
                    statement = line
                    line = ''
                else:
                    statement = line[:ind_semi]
                    line = line[ind_semi+1:]
                if statement: result = ((result + '\n') if result else '') + max_to_string(maxima_eval("#$%s$"%statement))
            else:
                statement = line[:ind_dollar]
                line = line[ind_dollar+1:]
                if statement: _ = maxima_eval("#$%s$"%statement)
        if not reformat:
            return result
        return ''.join([x.strip() for x in result.split()])

    _eval_line = eval

    ###########################################
    # Direct access to underlying lisp interpreter.
    ###########################################
    def lisp(self, cmd):
        """
        Send a lisp command to maxima.

        .. note::

           The output of this command is very raw - not pretty.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.lisp("(+ 2 17)")
            <ECL: 19>
        """
        return ecl_eval(cmd)

    def set(self, var, value):
        """
        Set the variable var to the given value.

        INPUT:


        -  ``var`` - string

        -  ``value`` - string


        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
        """
        if not isinstance(value, str):
            raise TypeError
        cmd = '%s : %s$'%(var, value.rstrip(';'))
        self.eval(cmd)

    def clear(self, var):
        """
        Clear the variable named var.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
            sage: maxima_lib.clear('xxxxx')
            sage: maxima_lib.get('xxxxx')
            'xxxxx'
        """
        try:
            self.eval('kill(%s)$'%var)
        except (TypeError, AttributeError):
            pass

    def get(self, var):
        """
        Get the string value of the variable var.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib.set('xxxxx', '2')
            sage: maxima_lib.get('xxxxx')
            '2'
        """
        s = self.eval('%s;'%var)
        return s

    def _create(self, value, name=None):
        name = self._next_var_name() if name is None else name
        if isinstance(value,EclObject):
            maxima_eval([[msetq],cadadr("#$%s$#$"%name),value])
        else:
            self.set(name, value)
        return name

    def _function_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._function_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibFunction'>
        """
        return MaximaLibFunction

    def _object_class(self):
        """
        Return the Python class of Maxima elements.

        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._object_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibElement'>
        """
        return MaximaLibElement

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._function_element_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibFunctionElement'>
        """
        return MaximaLibFunctionElement

    def _object_function_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: maxima_lib._object_function_class()
            <class 'sage.interfaces.maxima_lib.MaximaLibElementFunction'>
        """
        return MaximaLibElementFunction

    ##some helper functions to wrap the calculus use of the maxima interface.
    ##these routines expect arguments living in the symbolic ring and return something
    ##that is hopefully coercible into the symbolic ring again.

    def sr_integral(self,*args):
        try:
            return max_to_sr(maxima_eval(([max_integrate],[sr_to_max(SR(a)) for a in args])))
        except RuntimeError, error:
            s = str(error)
            if "Divergent" in s or "divergent" in s:
                raise ValueError, "Integral is divergent."
            else:
                raise error

    def sr_sum(self,*args):
        try:
            return max_to_sr(maxima_eval([[max_ratsimp],[[max_simplify_sum],([max_sum],[sr_to_max(SR(a)) for a in args])]]));
        except RuntimeError, error:
            s = str(error)
            if "divergent" in s:
                raise ValueError, "Sum is divergent."
            else:
                raise error

    def sr_limit(self,expr,v,a,dir=None):
        L=[sr_to_max(SR(a)) for a in [expr,v,a]]
        if dir == "plus":
            L.append(max_plus)
        elif dir == "minus":
            L.append(max_minus)
        return max_to_sr(maxima_eval(([max_limit],L)))

    def sr_tlimit(self,expr,v,a,dir=None):
        L=[sr_to_max(SR(a)) for a in [expr,v,a]]
        if dir == "plus":
            L.append(max_plus)
        elif dir == "minus":
            L.append(max_minus)
        return max_to_sr(maxima_eval(([max_tlimit],L)))


def is_MaximaLibElement(x):
    """
    Returns True if x is of type MaximaLibElement.

    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import maxima_lib, is_MaximaLibElement
        sage: m = maxima_lib(1)
        sage: is_MaximaLibElement(m)
        True
        sage: is_MaximaLibElement(1)
        False
    """
    return isinstance(x, MaximaLibElement)

class MaximaLibElement(MaximaAbstractElement):
    """
    """
    def ecl(self):
        try:
            return self._ecl
        except AttributeError:
            self._ecl=maxima_eval("#$%s$"%self._name)
            return self._ecl

    def to_poly_solve(self,vars,options=""):
        if options.find("use_grobner=true") != -1:
            cmd=EclObject([[max_to_poly_solve], self.ecl(), sr_to_max(vars),
                                             [[mequal],max_use_grobner,True]])
        else:
            cmd=EclObject([[max_to_poly_solve], self.ecl(), sr_to_max(vars)])
        return self.parent()(maxima_eval(cmd))

    def display2d(self, onscreen=True):
        """
        EXAMPLES::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: F = maxima_lib('x^5 - y^5').factor()
            sage: F.display2d ()
                                   4      3    2  2    3      4
                       - (y - x) (y  + x y  + x  y  + x  y + x )
        """
        self._check_valid()
        P = self.parent()
        P._eval_line('display2d : true$')
        s = stdout_to_string('disp(%s)'%self.name())
        #s = P._eval_line('disp(%s)$'%self.name())
        P._eval_line('display2d : false$')
        s = s.strip('\r\n')

        # if ever want to dedent, see
        # http://mail.python.org/pipermail/python-list/2006-December/420033.html
        if onscreen:
            print s
        else:
            return s


class MaximaLibFunctionElement(MaximaAbstractFunctionElement):
    pass


class MaximaLibFunction(MaximaAbstractFunction):
    pass


class MaximaLibElementFunction(MaximaLibElement, MaximaAbstractElementFunction):
    pass


# The (unique) instance
maxima_lib = MaximaLib()
maxima = maxima_lib


def reduce_load_MaximaLib():
    """
    EXAMPLES::

        sage: from sage.interfaces.maxima_lib import reduce_load_MaximaLib
        sage: reduce_load_MaximaLib()
        Maxima_lib
    """
    return maxima_lib

#**********************************
# ???

import sage.symbolic.expression
import sage.functions.trig
import sage.functions.log
import sage.functions.other
import sage.symbolic.integration.integral
from sage.symbolic.operators import FDerivativeOperator

car=EclObject("car")
cdr=EclObject("cdr")
caar=EclObject("caar")
cadr=EclObject("cadr")
cddr=EclObject("cddr")
caddr=EclObject("caddr")
caaadr=EclObject("caaadr")
cadadr=EclObject("cadadr")
meval=EclObject("meval")
NIL=EclObject("NIL")
ratdisrep=EclObject("ratdisrep")

sage_op_dict = {
    sage.symbolic.expression.operator.abs : "MABS",
    sage.symbolic.expression.operator.add : "MPLUS",
    sage.symbolic.expression.operator.div : "MQUOTIENT",
    sage.symbolic.expression.operator.eq : "MEQUAL",
    sage.symbolic.expression.operator.ge : "MGEQP",
    sage.symbolic.expression.operator.gt : "MGREATERP",
    sage.symbolic.expression.operator.le : "MLEQP",
    sage.symbolic.expression.operator.lt : "MLESSP",
    sage.symbolic.expression.operator.mul : "MTIMES",
    sage.symbolic.expression.operator.ne : "MNOTEQUAL",
    sage.symbolic.expression.operator.neg : "MMINUS",
    sage.symbolic.expression.operator.pow : "MEXPT",
    sage.symbolic.expression.operator.or_ : "MOR",
    sage.symbolic.expression.operator.and_ : "MAND",
    sage.functions.trig.acos : "%ACOS",
    sage.functions.trig.acot : "%ACOT",
    sage.functions.trig.acsc : "%ACSC",
    sage.functions.trig.asec : "%ASEC",
    sage.functions.trig.asin : "%ASIN",
    sage.functions.trig.atan : "%ATAN",
    sage.functions.trig.cos : "%COS",
    sage.functions.trig.cot : "%COT",
    sage.functions.trig.csc : "%CSC",
    sage.functions.trig.sec : "%SEC",
    sage.functions.trig.sin : "%SIN",
    sage.functions.trig.tan : "%TAN",
    sage.functions.log.exp : "%EXP",
    sage.functions.log.ln : "%LOG",
    sage.functions.log.log : "%LOG",
    sage.functions.other.factorial : "MFACTORIAL",
    sage.functions.other.erf : "%ERF",
    sage.functions.other.gamma_inc : "%GAMMA_INCOMPLETE"
}
#we compile the dictionary
sage_op_dict = dict([(k,EclObject(sage_op_dict[k])) for k in sage_op_dict])
max_op_dict = dict([(sage_op_dict[k],k) for k in sage_op_dict])

def add_vararg(*args):
    S=0
    for a in args:
        S=S+a
    return S

def mul_vararg(*args):
    P=1
    for a in args:
        P=P*a
    return P

def sage_rat(x,y):
    return x/y

mplus=EclObject("MPLUS")
mtimes=EclObject("MTIMES")
mdiff=EclObject("%DERIVATIVE")
rat=EclObject("RAT")
max_i=EclObject("$%I")
max_op_dict[mplus]=add_vararg
max_op_dict[mtimes]=mul_vararg
max_op_dict[rat]=sage_rat
mqapply=EclObject("MQAPPLY")
max_li=EclObject("$LI")
max_psi=EclObject("$PSI")
max_array=EclObject("ARRAY")
max_gamma_incomplete=sage_op_dict[sage.functions.other.gamma_inc]

def mrat_to_sage(expr):
    r"""
    Convert a maxima MRAT expression to Sage SR

    Maxima has an optimised representation for multivariate rational expressions.
    The easiest way to translate those to SR is by first asking maxima to give
    the generic representation of the object. That is what RATDISREP does in
    maxima.
    """
    return max_to_sr(meval(EclObject([[ratdisrep],expr])))

def mqapply_to_sage(expr):
    r"""
    Special conversion rule for MQAPPLY expressions
    """
    if caaadr(expr) == max_li:
        return sage.functions.log.polylog(max_to_sr(cadadr(expr)),max_to_sr(caddr(expr)))
    if caaadr(expr) == max_psi:
        return sage.functions.other.psi(max_to_sr(cadadr(expr)),max_to_sr(caddr(expr)))
    else:
        op=max_to_sr(cadr(expr))
        max_args=cddr(expr)
        args=[max_to_sr(a) for a in max_args]
        return op(*args)

def dummy_integrate(expr):
    r"""
    we would like to simply tie maxima's integrate to sage.calculus.calculus.dummy_integrate, but we're being imported there so to avoid circularity we define it here.
    """
    args=[max_to_sr(a) for a in cdr(expr)]
    if len(args) == 4 :
        return sage.symbolic.integration.integral.definite_integral(*args, hold=True)
    else:
        return sage.symbolic.integration.integral.indefinite_integral(*args, hold=True)

def mdiff_to_sage(expr):
    return max_to_sr(expr.cadr()).diff(*[max_to_sr(e) for e in expr.cddr()])

special_max_to_sage={
    EclObject("MRAT") : mrat_to_sage,
    mqapply : mqapply_to_sage,
    EclObject("%INTEGRATE") : dummy_integrate,
    mdiff : mdiff_to_sage
}

special_sage_to_max={
    sage.functions.log.polylog : lambda N,X : [[mqapply],[[max_li, max_array],N],X],
    sage.functions.other.psi1 : lambda X : [[mqapply],[[max_psi, max_array],0],X],
    sage.functions.other.psi2 : lambda N,X : [[mqapply],[[max_psi, max_array],N],X],
    sage.functions.other.Ei : lambda X : [[max_gamma_incomplete], 0, X]
}

sage_sym_dict={}
max_sym_dict={}

def pyobject_to_max(obj):
    if isinstance(obj,sage.rings.rational.Rational):
        return EclObject(obj) if (obj.denom().is_one()) else EclObject([[rat], obj.numer(),obj.denom()])
    elif isinstance(obj,sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic) and obj.parent().defining_polynomial().list() == [1,0,1]:
        re, im = obj.list()
        return EclObject([[mplus], pyobject_to_max(re), [[mtimes], pyobject_to_max(im), max_i]])
    return EclObject(obj)

# This goes from SR to EclObject
def sr_to_max(expr):
    r"""
    """
    global sage_op_dict, max_op_dict
    global sage_sym_dict, max_sym_dict
    if isinstance(expr,list) or isinstance(expr,tuple):
        return EclObject(([mlist],[sr_to_max(e) for e in expr]))
    op = expr.operator()
    if op:
        # Stolen from sage.symbolic.expression_conversion
        # Should be defined in a function and then put in special_sage_to_max
        # For that, we should change the API of the functions there (we need to have access to op, not only to expr.operands()
        if isinstance(op, FDerivativeOperator):
            from sage.symbolic.ring import is_SymbolicVariable
            args = expr.operands()
            if (not all(is_SymbolicVariable(v) for v in args) or
                len(args) != len(set(args))):
                raise NotImplementedError, "arguments must be distinct variables"
            f = sr_to_max(op.function()(*args))
            params = op.parameter_set()
            deriv_max = []
            [deriv_max.extend([sr_to_max(args[i]), EclObject(params.count(i))]) for i in set(params)]
            l = [mdiff,f]
            l.extend(deriv_max)
            return EclObject(l)
        elif (op in special_sage_to_max):
            return EclObject(special_sage_to_max[op](*[sr_to_max(o) for o in expr.operands()]))
        elif not (op in sage_op_dict):
            # Maxima does some simplifications automatically by default
            # so calling maxima(expr) can change the structure of expr
            #op_max=caar(maxima(expr).ecl())
            # This should be safe if we treated all special operators above
            op_max=maxima(op).ecl()
            sage_op_dict[op]=op_max
            max_op_dict[op_max]=op
        return EclObject(([sage_op_dict[op]],
                     [sr_to_max(o) for o in expr.operands()]))
    elif expr._is_symbol() or expr._is_constant():
        if not expr in sage_sym_dict:
            sym_max=maxima(expr).ecl()
            sage_sym_dict[expr]=sym_max
            max_sym_dict[sym_max]=expr
        return sage_sym_dict[expr]
    else:
        try:
            return pyobject_to_max(expr.pyobject())
        except TypeError:
            return maxima(expr).ecl()

# This goes from EclObject to SR
def max_to_sr(expr):
    if expr.consp():
        op_max=caar(expr)
        if op_max in special_max_to_sage:
            return special_max_to_sage[op_max](expr)
        if not(op_max in max_op_dict):
            # This could be unsafe if the conversion to SR chenges the structure of expr
            sage_expr=SR(maxima(expr))
            max_op_dict[op_max]=sage_expr.operator()
            sage_op_dict[sage_expr.operator()]=op_max
        op=max_op_dict[op_max]
        max_args=cdr(expr)
        args=[max_to_sr(a) for a in max_args]
        return op(*args)
    elif expr.symbolp():
        if not(expr in max_sym_dict):
            sage_symbol=SR(maxima(expr))
            sage_sym_dict[sage_symbol]=expr
            max_sym_dict[expr]=sage_symbol
        return max_sym_dict[expr]
    else:
        return expr.python()

