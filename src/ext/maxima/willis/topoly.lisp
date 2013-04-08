;; Written by Barton Willis, willisb@unk.edu
;; and emailed to maxima users group email list
;; on Sep 2006
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defun suppress-multiple-zeros (q)
  (let ((acc 1) ($factorflag nil))
    (setq q ($factor q))
    (setq q (if (mtimesp q) (margs q) (list q)))
    (dolist (qi q acc)
      (setq acc (mul acc (cond ((mnump qi) (if (eq t (meqp qi 0)) 0 1))
			       ((mexptp qi) (nth 1 qi))
			       (t qi)))))))

;; When the second argument is true, convert constant expressions to polynomial form.
;; The default for the second argument is false. Examples:
;;
;;(%i1) topoly(x-sqrt(2),false);
;;(%o1) x-sqrt(2)=0
;;(%i2) topoly(x-sqrt(2),true);
;;(%o2) x^2-2=0
;;(%i7) topoly(sqrt(2)+sqrt(3)-x,true);
;;(%o7) x^4-10*x^2+1=0
;;(%i8) solve(%,x);
;;(%o8) [x=-sqrt(2*sqrt(6)+5),x=sqrt(2*sqrt(6)+5),x=-sqrt(5-2*sqrt(6)),x=sqrt(5-2*sqrt(6))]
;;(%i9) map(lambda([e],topoly(e,true)),%);
;;(%o9) [x^4-10*x^2+1=0,x^4-10*x^2+1=0,x^4-10*x^2+1=0,x^4-10*x^2+1=0]

(defun $topoly (p &optional (cnst nil))

  (if (not (or (eq t cnst) (eq nil cnst)))
      (merror "The second argument to 'topoly' must be either 'true' or 'false'"))

  (let ((subs) (q) (nv `((mlist)))) ;; new variables
    (setq p (meqhk p))
    (setq q ($ratdenom p))
    (if (not ($constantp q)) (mtell "Assuming that ~:M " `((mnotequal) ,q 0)))
    (setq p ($ratdisrep ($ratnumer p)))

    (setq p (to-polynomial p nil cnst))
    (setq subs (second p))
    (setq p (first p))
    (dolist (sk subs)
      (setq nv ($append nv ($listofvars ($lhs sk)))))
    (setq p (if (null subs) p ($first (mfuncall '$eliminate `((mlist) ,p ,@subs) nv))))
    `((mequal) ,(suppress-multiple-zeros p) 0)))

(defun to-polynomial (p subs cnst)
  (cond (($mapatom p) (list p subs))
	((and (not cnst) ($constantp p)) (list p subs))

	((mexptp p)
	 (let ((n (nth 2 p)) (b (nth 1 p)) (nv)(l))
	   (cond ((integerp n)
		  (setq b (to-polynomial b nil cnst))
		  (setq subs (append (second b) subs))
		  (setq b (first b))
		  (if (> n 0) (list (power b n) subs) (merror "Unable to convert to a polynomial equation")))

		 (($ratnump n)

		  (setq b (to-polynomial b nil cnst))
		  (setq subs (append (second b) subs))
		  (setq b (first b))
		  (setq nv (gensym))
		  (setq subs (cons `((mequal) ,(power nv ($denom n)) ,(power b ($num n))) subs))
		  (list nv subs))
		 (t (merror "Nonalgebraic argument given to 'topoly'")))))

	((op-equalp p 'mabs)
	 (setq b (to-polynomial (first (margs p)) nil cnst))
	 (setq subs (append (second b) subs))
	 (setq b (first b))
	 (setq nv (gensym))
	 (list nv (cons `((mequal) ,(power nv 2) ,(power b 2)) subs)))

	((mtimesp p)
	 (let ((z 1) (acc nil))
	   (setq p (mapcar #'(lambda (s) (to-polynomial s nil cnst)) (margs p)))
	   (dolist (pk p)
	     (setq z (mul z (first pk)))
	     (setq acc (append acc (second pk))))
	   (list z acc)))

	 ((mplusp p)
	  (let ((z 0) (acc nil))
	    (setq p (mapcar #'(lambda (s) (to-polynomial s nil cnst)) (margs p)))
	    (dolist (pk p)
	      (setq z (add z (first pk)))
	      (setq acc (append acc (second pk))))
	    (list z acc)))

	 (t (merror "Nonalgebraic argument given to 'topoly'"))))





