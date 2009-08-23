;;  Author Barton Willis
;;  University of Nebraska at Kearney
;;  Copyright (C) 2006, 2007, 2008 Barton Willis

;;  This program is free software; you can redistribute it and/or modify
;;  it under the terms of the GNU General Public License as published by
;;  the Free Software Foundation; either version 2 of the License, or
;;  (at your option) any later version.

;;  This program is distributed in the hope that it will be useful,
;;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;  GNU General Public License for more details.


($put '$to_poly 2 '$version)
($load '$polynomialp)

(defmacro opapply (op args)
  `(simplify (cons (list ,op) ,args)))

;; The next three functions convert max and min to abs functions.

(defun max-to-abs (e)
  (reduce #'(lambda (a b) (div (add (add a b) (take '(mabs) (sub a b))) 2)) e))

(defun min-to-abs (e)
  (reduce #'(lambda (a b) (div (sub (add a b) (take '(mabs) (sub a b))) 2)) e))

(defun convert-from-max-min-to-abs (e)
  (cond (($mapatom e) e)
	((op-equalp e '$max) (max-to-abs (mapcar 'convert-from-max-min-to-abs (margs e))))
	((op-equalp e '$min) (min-to-abs (mapcar 'convert-from-max-min-to-abs (margs e))))
	(t (opapply (mop e) (mapcar 'convert-from-max-min-to-abs (margs e))))))

(defun maxima-variable-p (e)
  (or (symbolp e) ($subvarp e)))

(defun list-subst (l p)
  (if (null l) p (list-subst (rest l) ($substitute (first l) p))))

;; to_poly(p,vars) returns a polynomial in the variables 'vars' that has a zero whenever
;; p has a zero. When 1 is a member of vars, constant terms, such as sqrt(5) also get
;; converted to polynomial form. The value of vars defaults to all variables including
;; constants.

(defun $to_poly (p &optional (vars 'convert-all-vars))
  (let (($listconstvars t) (q) (convert-cnst nil) (proviso nil) (non-alg) (g-vars nil))

    (if (eq vars 'convert-all-vars) (setq vars ($cons 1 ($listofvars p))))

    (if (not ($listp vars))
	(merror "The second argument to 'topoly' must be a list"))

    (cond (($member 1 vars)
	   (setq convert-cnst t)
	   (setq vars ($delete 1 vars))))
    (setq p (meqhk p))
    (setq q ($ratdenom p))
    (if (not ($constantp q)) (push `((mnotequal) ,(sratsimp q) 0) proviso))
    (multiple-value-setq (p g-vars) (non-algebraic-subst (sratsimp ($ratnumer p)) vars))
    (setq non-alg ($second p))
    (setq p ($first p))

    (setq q ($ratdenom p))
    (if (not ($constantp q)) (push `((mnotequal) ,(sratsimp q) 0) proviso))
    (setq p (sratsimp ($ratnumer p)))

    (push '(mlist) g-vars)
    (setq vars ($append vars g-vars))

    ;; It's OK to send every expression through convert-from-max-min-to-abs.
    ;; I put in the conditional to skip the ratsimp for expressions that don't
    ;; involve max or min.

    (setq p (if ($freeof '$max '$min p) p (sratsimp (convert-from-max-min-to-abs p))))
    (setq p (to-polynomial p vars convert-cnst))
    (setq proviso (opapply 'mlist (append proviso (third p))))
    `((mlist) ((mlist) ,(first p) ,@(second p)) ,proviso ,non-alg)))

(defun to-polynomial (p vars convert-cnst)
  (let ((n) (b) (nv) (acc nil) (subs nil) (pk) (q) (inequal) (np-subs))
    (cond ((or (maxima-variable-p p)
	       (mnump p)
	       (and ($emptyp vars) (not convert-cnst))
	       (and (not ($constantp p)) ($lfreeof vars p))
	       (and ($constantp p) (not convert-cnst))
	       (and (op-equalp p '$conjugate) (maxima-variable-p (second p))))
	   (list p nil nil nil))

	  ((mexptp p)

	   (setq n (nth 2 p))
	   (setq b (nth 1 p))
	   (cond ((and (integerp n) (> n 0))
		  (list p nil nil nil))

		 (($ratnump n)
		  (setq b (to-polynomial b vars convert-cnst))
		  (setq subs (second b))
		  (setq inequal (third b))
		  (setq np-subs (fourth b))
		  (setq b (first b))
		  (setq nv (new-gentemp 'general))
		  (cond ((or (mgrp n 0) (mnump b))
			 (let ((q) (r) (rr))
			   (setq q (take '($floor) n))
			   (setq r (sub n q))
			   (setq rr ($denom r))
			   (push (take '(mleqp) (take '($parg) nv) (div '$%pi rr)) inequal)
			   (push (take '(mlessp) (div '$%pi (neg rr)) (take '($parg) nv)) inequal)
			   (push (take '(mequal) (power b ($num r)) (power nv ($denom r))) subs)
			   (push (take '(mequal) p nv) np-subs)
			   (list (mul nv (power b q))  subs inequal np-subs)))

			;   (push (take '(mleqp) (take '($parg) nv) (mul n '$%pi)) inequal)
			;   (push (take '(mlessp) (mul -1 '$%pi n) (take '($parg) nv)) inequal)
			;   (push (take '(mequal) (power b ($num n)) (power nv ($denom n))) subs)
			;   (push (take '(mequal) p nv) np-subs)
			;   (list nv subs inequal np-subs)))

			(t
			 (setq n (neg n))
			 (push (take '(mequal) 1 (mul (power nv ($denom n)) (power b ($num n)))) subs)
			 (push (take '(mlessp) (take '($parg) nv) (mul n '$%pi)) inequal)
			 (push (take '(mleqp) (mul -1 '$%pi n) (take '($parg) nv)) inequal)
			 (push (take '(mnotequal) nv 0) inequal)
			 (list nv subs inequal np-subs))))

		 (t (merror "Nonalgebraic argument given to 'topoly'"))))

	  ((op-equalp p 'mabs)
	   (setq b (to-polynomial (first (margs p)) vars convert-cnst))
	   (setq acc (second b))
	   (setq inequal (third b))
	   (setq np-subs (fourth b))
	   (setq b (first b))
	   (setq nv (new-gentemp '$general))
	   ;; Ok, I'm indecisive.
	   ;(push (take '(mgeqp) nv 0) inequal)
	   ;(push (take '(mequal) (take '($parg) nv) 0) inequal)
	   (push (take '($isnonnegative_p) nv) inequal)
	   (list nv (cons (take '(mequal) (mul b (take '($conjugate) b)) (mul nv (take '($conjugate) nv))) acc)
		 inequal np-subs))

	  ((mtimesp p)
	   (setq acc 1)
	   (setq p (margs p))
	   (while p
	     (setq pk (first p))
	     (setq p (rest p))
	     (setq q (to-polynomial pk vars convert-cnst))
	     (setq acc (mul acc (first q)))
	     (setq subs (append (second q) subs))
	     (setq inequal (append inequal (third q)))
	     (setq np-subs (append np-subs (fourth q)))
	     (setq vars ($append vars ($listofvars `((mlist) ,@subs))))

	     (setq p (mapcar #'(lambda (s) (list-subst np-subs s)) p)))
	   (list acc subs inequal np-subs))

	  ((mplusp p)
	   (setq acc 0)
	   (setq p (margs p))
	   (while p
	     (setq pk (first p))
	     (setq p (rest p))
	     (setq q (to-polynomial pk vars convert-cnst))
	     (setq acc (add acc (first q)))
	     (setq subs (append (second q) subs))
	     (setq inequal (append (third q) inequal))
	     (setq np-subs (append (fourth q) np-subs))
	     (setq vars ($append vars ($listofvars `((mlist) ,@subs))))
	     (setq p (mapcar #'(lambda (s) (list-subst np-subs s)) p)))
	   (list acc subs inequal np-subs))


	  (t (merror "Nonalgebraic argument given to 'topoly'")))))


#|
  Things I don't like about eliminate:

(1)  eliminate([x + x*y + z-4, x+y+z-12, x^2 + y^2 + z^2=7],[x,y,z]) -> [4]

Here Maxima solves for z. There are more than one solution for z, but Maxima returns
just one solution. A user might think that there is one solution or that
the equations are inconsistent.

(2)  eliminate([x+y-1,x+y-8],[x,y]) -> Argument to `last' is empty.

Here the equations are inconsistent, but we get an error (from solve) instead
of a clear message about what happened.

(3) eliminate([a],[x]) -> Can't eliminate from only one equation -- an error.
but eliminate([a,a],[x]) -> [a,a]. This is silly.

(4) eliminate([x],[]) -> Can't eliminate from only one equation.
but eliminate([x,x],[]) -> [x,x]. Again, this is silly.

(5) elim([x^3-y^2,x^7 + y+z*p,q*q-23+ x+y,p-q,x+y+z-5],[x,y,z,p]) takes 0.3 second
but eliminate([x^3-y^2,x^7 + y+z*p,q*q-23+ x+y,p-q,x+y+z-5],[x,y,z,p]) takes
a long time (I never let it run to completion). Unlike 'eliminate,' the function 'elim'
makes some guesses about which polynomial to use as the pivot and which variable
to eliminate.
|#

(defun require-maxima-variable (x context-string)
  (setq x (ratdisrep x))
  (if (maxima-variable-p x) x
    (merror "Function ~:M expects a symbol, instead found ~:M" context-string x)))

;; Simplify a polynomial equation p = 0 by

;;   (1) nonzero constant * p --> p,
;;   (2) p^n --> p.

;; If you want to use $factor instead of $sqfr, go ahead. But if you do that, you might want to
;; locally set factorflag to false.

(defun suppress-multiple-zeros (q)
  (setq q ($sqfr q))
  (cond ((mnump q) (if (zerop1 q) 0 1))
	((mtimesp q) (muln (mapcar 'suppress-multiple-zeros (margs q)) t))
	((and (mexptp q) (integerp (third q)) (> (third q) 0)) (second q))
	(($constantp q) 1) ; takes care of things like (1 + sqrt(5)) * x --> x.
	(t q)))

;; Using eq as the "pivot," eliminate x from the list or set of equations eqs.

(defun $eliminate_using (eqs eq x)
  (if (or ($listp eqs) ($setp eqs))
      (progn
	(setq eqs (mapcar #'(lambda (s) ($ratexpand (meqhk s))) (margs eqs)))
	(setq eqs (margs (opapply '$set eqs))))
    (merror "The first argument to 'eliminate_using' must be a list or a set"))

  (setq x (require-maxima-variable x "$eliminate_using"))

  (if (not (every #'(lambda (s) ($polynomialp s `((mlist) ,x)
					      `((lambda) ((mlist) s) (($freeof) ,x s)))) eqs))
      (merror "The first argument to 'eliminate_using' must be a set or list of polynomials"))

  (setq eq ($ratexpand (meqhk eq)))
  (if (not ($polynomialp eq `((mlist) ,x) `((lambda) ((mlist) s) (($freeof) ,x s))))
      (merror "The second argument to 'eliminate_using' must be a polynomial"))

  (setq eqs (mapcar #'suppress-multiple-zeros eqs))
  (setq eq (suppress-multiple-zeros eq))
  (opapply '$set (eliminate-using eqs eq x)))

(defun eliminate-using (eqs eq x)
  (delete 0 (mapcar #'(lambda (s) (suppress-multiple-zeros ($resultant s eq x))) eqs)))

;; Return an upper bound for the total degree of the polynomial p in the variables vars.
;; When p is fully expanded, the bound is exact.

(defun degree-upper-bound (p vars)
  (cond ((maxima-variable-p p) (if (member p vars :test #'equal) 1 0))

	((mnump p) 0)

	((and (mexptp p) (integerp (third p)) (> (third p) 0))
	 (* (degree-upper-bound (nth 1 p) vars) (nth 2 p)))

	((mtimesp p)
	 (reduce '+ (mapcar #'(lambda (s) (degree-upper-bound s vars)) (margs p))))

	((mplusp p)
	 (simplify `(($max) ,@(mapcar #'(lambda (s) (degree-upper-bound s vars)) (margs p)))))

	((apply '$freeof (append vars (list p))) 0)
	(t (merror "Nonpolynomial argument given to degree-upper-bound"))))

(defun unk-eliminate (eqs vars &optional (pivots nil))
   (let ((ni) (n-min nil) (xeqs `(($set))) (pivot-var) (pivot-eq) (acc `(($set))) ($ratfac nil))
     (cond ((or (null vars) (null eqs)) (list eqs pivots))
	   (t

	    ;; The pivot is a nonconstant member of eqs with minimal total degree.
	    ;; The  constant members of eqs get adjoined into acc -- all other  members get
	    ;; adjoined into xeqs. Since each member of eqs has been ratexpanded,
	    ;; the degree bound is exact.

	    (dolist (ei eqs)
	      (setq ei ($ratexpand (suppress-multiple-zeros ei)))
	      (setq ni (degree-upper-bound ei vars))

	      (if (and (or (null n-min) (< ni n-min)) (> ni 0))
		  (setq n-min ni pivot-eq ei))

	      (if (and (> ni 0) (not (equal 0 ei))) (setq xeqs ($adjoin ei xeqs)) (setq acc ($adjoin ei acc))))

	    (setq xeqs (margs xeqs))
	    (setq acc (margs acc))
	    ;; Now we'll decide which variable to eliminate. The pivot variable
	    ;; is the variable that has the least (but nonzero) degree in pivot-eq.

	    (setq n-min nil)
	    (dolist (vi vars)
	      (setq ni (degree-upper-bound pivot-eq (list vi)))
	      (if (and (or (null n-min) (< ni n-min)) (> ni 0))
		  (setq pivot-var vi n-min ni)))

	    (if (null pivot-var) (list eqs pivots)
	      (unk-eliminate (append acc (eliminate-using xeqs pivot-eq pivot-var)) (delete pivot-var vars)
			     (cons pivot-eq pivots)))))))

(defun $elim (eqs x)
  (if (or ($listp eqs) ($setp eqs))
      (progn
	(setq eqs (mapcar #'(lambda (s) ($ratexpand (suppress-multiple-zeros (meqhk s)))) (margs eqs)))
	(setq eqs (margs (opapply '$set eqs))))
    (merror "The first argument to 'elim' must be a list or a set"))

  (setq x (margs (cond (($listp x) ($setify x))
		       (($setp x) x)
		       (t (merror "The second argument to 'elim' must be a list or a set")))))

  (setq x (mapcar #'(lambda (s) (require-maxima-variable s "$elim")) x))

  (setq x (opapply 'mlist x))
  (if (not (every #'(lambda (s) ($polynomialp s x `((lambda) ((mlist) s) (($lfreeof) ,x s)))) eqs))
      (merror "Each member of the first argument to 'elim' must be a polynomial"))

  (setq x (margs x))
  (opapply 'mlist (mapcar #'(lambda (s) (opapply 'mlist s)) (unk-eliminate eqs x))))

(defun $elim_allbut (eqs x)
  (let (($listconstvars nil) (v))
    (setq v ($listofvars eqs))
    (setq x (cond (($listp x) ($setify x))
		  (($setp x) x)
		  (t (take '($set) x))))
    (setq v ($setdifference ($setify v) x))
    ($elim eqs ($listify v))))

;; Return a list of the arguments to the operator op in the expression e.
;; This function only looks inside + and * and it doesn't recurse inside
;; + and *. So

;;  gather-args(log(x) + log(a + log(b)), log) --> (x, a + log(b)),
;;  gather-args(cos(log(x)), log) --> (),
;;  gather-args(x^log(x) + log(s)^x, log) --> ().

(defun gather-args (e op vars)
  (cond (($mapatom e) nil)
	((and (eq (mop e) op) (not ($lfreeof vars e))) (margs e))
	((memq (mop e) '(mplus mtimes))
	 (let ((acc nil))
	   (setq e (margs e))
	   (dolist (ek e acc)
	     (setq acc (append acc (gather-args ek op vars))))))
	(t nil)))

(defun gather-exp-args (e vars)
  (cond (($mapatom e) nil)
	((eq (mop e) 'mexpt)
	 (if (and (eq (second e) '$%e) (not ($lfreeof vars (third e)))) (list (third e))
	   (gather-exp-args (second e) vars)))
	(t (mapcan #'(lambda (s) (gather-exp-args s vars)) (margs e)))))

;; Return true iff the set {a,b} is linearly dependent. Thus return true iff
;; either a = 0 or b = 0, or a/b is a number.

(defun $linearly_dependent_p (a b)
  (linearly-dependent-p a b))

(defun linearly-dependent-p (a b)
  (if (or (=0 a)
	  (=0 b)
	  (mnump (sratsimp (div a b)))) t nil))

(defun exponentialize-if (e vars)
  (cond (($mapatom e) e)
	((and (trigp (mop e)) (not ($lfreeof vars e))) ($exponentialize e))
	(t (opapply (mop e) (mapcar #'(lambda (s) (exponentialize-if s vars)) (margs e))))))

(defun logarc-if (e vars)
  (cond (($mapatom e) e)
	((and (arcp (mop e)) (not ($lfreeof vars e))) ($logarc e))
	(t (opapply (mop e) (mapcar #'(lambda (s) (logarc-if s vars)) (margs e))))))

(defun non-algebraic-subst (e vars)
  (let ((log-args nil) (exp-args nil) (s) (g) (p) (q) (na-subs nil) (g-vars nil) (sz))

    (setq e ($ratexpand (exponentialize-if (logarc-if e vars) vars)))
    (setq log-args (gather-args e '%log vars))
    (setq log-args (margs (opapply '$set log-args)))
    (setq exp-args (gather-exp-args e vars))
    (setq exp-args (opapply '$set exp-args))
    (setq s (margs (simplify ($equiv_classes exp-args 'linearly-dependent-p))))

    (dolist (sk s)
      (setq g (new-gentemp 'general))
      (push g g-vars)
      (setq sz ($first (apply '$ezgcd (margs sk))))
      (setq p (power '$%e sz))
      (push (take '(mequal) g p) na-subs)
      (setq sk (margs sk))
      (dolist (sl sk)
	(setq q (sratsimp (div sl sz)))
	(setq e ($ratexpand ($substitute (power g  q) (power '$%e sl) e)))))

    (dolist (sk log-args)
      (setq g (new-gentemp 'general))
      (push g g-vars)
      (setq p (take '(%log) sk))
      (setq e ($substitute g p e))
      (push (take '(mequal) g p) na-subs))

    (values `((mlist) ,e ((mlist) ,@na-subs)) g-vars)))

;; A simplifying carg function.

(setf (get '$parg 'operators) 'simp-parg)

(defun simp-parg (e yy z)
  (declare (ignore yy))
  (let (($domain '$complex) (dosimp nil) (complexsign t) (ec) (x) (y) (x-sgn) (y-sgn))
    (oneargcheck e)
    (setq e (simplifya (specrepcheck (cadr e)) z))

    (cond ((and (mtimesp e) ($numberp (cadr e))) ;; number * product
	   (cond ((eq t (mgrp (cadr e) 0))
		  (take '($parg) (reduce 'mul (cddr e))))

		 ((eq t (mgrp 0 (cadr e)))
		  (setq x (add '$%pi (take '($parg) (reduce 'mul (cddr e)))))
		  (sub '$%pi (simplify `(($mod) ,(sub '$%pi x) ,(mul 2 '$%pi)))))

		 (t `(($parg simp) ,e))))

	  ;; parg(constant^x) = parg(imagpart(x * log(constant)))
	  ((and (mexptp e) ($constantp (cadr e)) (eq t (mgrp (cadr e) 0)))
	   (setq x (mul (take '(%log) (cadr e)) (caddr e)))
	   (setq x (div (sub x (take '($conjugate) x)) (mul 2 '$%i)))
	   (sub '$%pi (simplify `(($mod) ,(sub '$%pi x) ,(mul 2 '$%pi)))))

	  ((and (mexptp e) ($constantp (caddr e)) (eq t (mgrp 0 (caddr e))) (eq t (mgrp 1 (caddr e))))
	   (mul (caddr e) (take '($parg) (cadr e))))

	  ((op-equalp e '%log)
	   (take '($parg) (cadr e)))

	  (t
	   (setq ec (take '($conjugate) e))
	   (setq x (div (add e ec) 2)) ;; the real part
	   (setq y (div (sub e ec) (mul 2 '$%i))) ;; the imaginary part

	   (setq x-sgn (csign x))
	   (cond ((eq x-sgn '$zero)
		  (mul (div '$%pi 2) (take '(%signum) y))) ;; imaginary axis (include pole).
		 (t
		  (setq y-sgn (csign y))
		  (cond ((eq x-sgn '$neg)
			 (cond ((eq y-sgn '$neg) ;; 3rd quadrant
				(sub (take '(%atan) (div y x)) '$%pi))
			       ((memq y-sgn '($zero $pz $pos)) ;; 2nd quadrant + negative real.
				(add (take '(%atan) (div y x)) '$%pi))
			       (t `(($parg simp) ,e))))

			((eq y-sgn '$zero) ;; real axis
			 (cond ((eq x-sgn '$neg) '$%pi)
			       ((memq x-sgn '($zero $pz $pos)) 0) ;; $zero shouldn't happen.
			       (t `(($parg simp) ,e))))

			((eq x-sgn '$pos) (take '(%atan) (div y x))) ;; 1st and 4th quadrants.
			(t  `(($parg simp) ,e)))))))))


(setf (get '$isreal_p 'operators) 'simp-isreal-p)

(defun simp-isreal-p (e yy z)
  (declare (ignore yy))
  (oneargcheck e)
  (let ((ec) ($domain '$complex))

    (setq e (simplifya (specrepcheck (cadr e)) z))

    ;; Simplifications:

    ;;  (1) if r is real then isreal_p(r + x) --> isreal_p(x),
    ;;  (2) if r is real then isreal_p(r * x) --> isreal_p(x),
    ;;  (3) if e = conjugate(e) then true,
    ;;  (4) if is(notequal(e,conjugate(e))) then false.

    (cond ((mtimesp e)
	   (setq e (muln (mapcar #'(lambda (s) (if (eq t (take '($isreal_p) s)) 1 s)) (margs e)) t)))
	  ((mplusp e)
	   (setq e (addn (mapcar #'(lambda (s) (if (eq t (take '($isreal_p) s)) 0 s)) (margs e)) t))))

    ;; If e is constant, apply rectform.
    (if ($constantp e) (setq e ($rectform e)))
    (setq ec (take '($conjugate) e))
    (cond ((eq t (meqp e ec)) t)
	  ((eq t (mnqp e ec)) nil)
	  (t `(($isreal_p simp) ,e)))))

(defvar *integer-gentemp-prefix* "$%Z")
(defvar *natural-gentemp-prefix* "$%N")
(defvar *real-gentemp-prefix*    "$%R")
(defvar *complex-gentemp-prefix* "$%C")
(defvar *general-gentemp-prefix* "$%G")

;; Return a new gentemp with a prefix that is determined by 'type.' Also put
;; an identifying tag on the symbol property of each new gentemp.

(defun $new_variable (type)
  (new-gentemp type))

(defun new-gentemp (type)
  (let ((g))
    (cond
     ((eq type '$integer)
      (setq g (gentemp *integer-gentemp-prefix*))
      (setf (get g 'integer-gentemp) t)
      (mfuncall '$declare g '$integer))

     ((eq type '$natural_number)
      (setq g (gentemp *natural-gentemp-prefix*))
      (setf (get g 'natural-gentemp) t)
      (mfuncall '$declare g '$integer)
      (mfuncall '$assume (take '(mgeqp) g 0)))

     ((eq type '$real)
      (setq g (gentemp *real-gentemp-prefix*))
      (setf (get g 'real-gentemp) t))

     ((eq type '$complex)
      (setq g (gentemp *complex-gentemp-prefix*))
      (setf (get g 'complex-gentemp) t)
      (mfuncall '$declare g '$complex))

     (t
      (setq g (gentemp *general-gentemp-prefix*))
      (setf (get g 'general-gentemp) t)))

    g))

;; Find all the new-gentemp variables in an expression e and re-index them starting from 0.

(defun $nicedummies (e)
  (let (($listconstvars nil)
	(z-vars nil) (z-k 0)
	(n-vars nil) (n-k 0)
	(r-vars nil) (r-k 0)
	(c-vars nil) (c-k 0))

    (mapcar #'(lambda (s) (let ((a))
			    (cond ((get s 'integer-gentemp)
				   (setq a (concat *integer-gentemp-prefix* z-k))
				   (mfuncall '$declare a '$integer)
				   (incf z-k)
				   (push `((mequal) ,s ,a) z-vars))

				  ((get s 'natural-gentemp)
				   (setq a (concat *natural-gentemp-prefix* n-k))
				   (mfuncall '$declare a '$integer)
				   (mfuncall '$assume (take '(mgeqp) a 0))
				   (incf n-k)
				   (push `((mequal) ,s ,a) n-vars))

				  ((get s 'real-gentemp)
				   (setq a (concat *real-gentemp-prefix* r-k))
				   (incf r-k)
				   (push `((mequal) ,s ,a) r-vars))

				  ((get s 'complex-gentemp)
				   (setq a (concat *complex-gentemp-prefix* c-k))
				   (mfuncall '$declare a '$complex)
				   (incf c-k)
				   (push `((mequal) ,s ,a) c-vars)))))
	    (margs ($sort ($listofvars e))))

    (push '(mlist) z-vars)
    (push '(mlist) n-vars)
    (push '(mlist) r-vars)
    (push '(mlist) c-vars)
    ($sublis c-vars ($sublis r-vars ($sublis n-vars ($sublis z-vars e))))))

(defun $complex_number_p (e)
  (complex-number-p e #'$numberp))
