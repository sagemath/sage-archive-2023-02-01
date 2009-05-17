;;;; `Pyrex' mode.

;; We define `pyrex-mode', a derived mode of `python-mode'.  We add some
;; font-lock keywords for things like cdef and cpdef, and we redefine some
;; `python-' functions to respect the additional pyrex keywords where
;; possible.  The advice is installed to not override in pure `python-mode'
;; buffers.

(require 'python)

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.pyx\\'" . pyrex-mode))
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.pxi\\'" . pyrex-mode))
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.pxd\\'" . pyrex-mode))

;;;###autoload
(define-derived-mode pyrex-mode python-mode "Pyrex"
  (set (make-local-variable 'outline-regexp)
       (rx (* space) (or "class" "def" "cdef" "cpdef" "elif" "else" "except" "finally"
			 "for" "if" "try" "while" "with")
	   symbol-end))

  (font-lock-add-keywords
   nil
   `((,(concat "\\<\\(NULL"
	       "\\|c\\(pdef\\|def\\|har\\|typedef\\|import\\)"
	       "\\|e\\(num\\|xtern\\)"
	       "\\|float"
	       "\\|in\\(clude\\|t\\)"
	       "\\|object\\|public\\|struct\\|type\\|union\\|void"
	       "\\)\\>")
      1 font-lock-keyword-face t))))

;; overload some python functions to better handle pyrex code
(defun pyrex-mode-p ()
  "Return t if we're in a pyrex-mode buffer."
  (derived-mode-p 'pyrex-mode))

(defvar python-font-lock-keywords
  `(,(rx symbol-start
	 ;; From v 2.5 reference, § keywords.
	 ;; def and class dealt with separately below
	 (or "and" "as" "assert" "break" "continue" "del" "elif" "else"
	     "except" "exec" "finally" "for" "from" "global" "if"
	     "import" "in" "is" "lambda" "not" "or" "pass" "print"
	     "raise" "return" "try" "while" "with" "yield"
             ;; Not real keywords, but close enough to be fontified as such
             "self" "True" "False")
	 symbol-end)
    (,(rx symbol-start "None" symbol-end) ; See § Keywords in 2.5 manual.
     . font-lock-constant-face)
    ;; Definitions
    (,(rx symbol-start (group "class") (1+ space) (group (1+ (or word ?_))))
     (1 font-lock-keyword-face) (2 font-lock-type-face))
    (,(rx symbol-start (group (or "def" "cdef" "cpdef")) (1+ space) (group (1+ (or word ?_))))
     (1 font-lock-keyword-face) (2 font-lock-function-name-face))
    ;; Top-level assignments are worth highlighting.
    (,(rx line-start (group (1+ (or word ?_))) (0+ space) "=")
     (1 font-lock-variable-name-face))
    (,(rx "@" (1+ (or word ?_))) ; decorators
    (0 font-lock-preprocessor-face))))

(defun pyrex-open-block-statement-p (&optional bos)
  "Return non-nil if statement at point opens a Pyrex block.
BOS non-nil means point is known to be at beginning of statement."
  (save-excursion
    (unless bos (python-beginning-of-statement))
    (looking-at (rx (and (or "if" "else" "elif" "while" "for" "def" "cdef" "cpdef"
			     "class" "try" "except" "finally" "with")
			 symbol-end)))))

(defadvice python-open-block-statement-p
  (around python-open-block-statement-p-advice first (&rest rest) activate)
  ;; this strange incantation calls the original python-open-block-statement
  ;; unless we're in a derived mode.  The (setq ad-return-value ...) is how
  ;; one modifies the return value of advised functions.
  (if (not (pyrex-mode-p))
      ad-do-it
    (setq ad-return-value (apply 'pyrex-open-block-statement-p rest))))

(defun pyrex-beginning-of-defun ()
  "`beginning-of-defun-function' for Pyrex.
Finds beginning of innermost nested class or method definition.
Returns the name of the definition found at the end, or nil if
reached start of buffer."
  (let ((ci (current-indentation))
	(def-re (rx line-start (0+ space) (or "def" "cdef" "cpdef" "class") (1+ space)
		    (group (1+ (or word (syntax symbol))))))
	found lep) ;; def-line
    (if (python-comment-line-p)
	(setq ci most-positive-fixnum))
    (while (and (not (bobp)) (not found))
      ;; Treat bol at beginning of function as outside function so
      ;; that successive C-M-a makes progress backwards.
      ;;(setq def-line (looking-at def-re))
      (unless (bolp) (end-of-line))
      (setq lep (line-end-position))
      (if (and (re-search-backward def-re nil 'move)
	       ;; Must be less indented or matching top level, or
	       ;; equally indented if we started on a definition line.
	       (let ((in (current-indentation)))
		 (or (and (zerop ci) (zerop in))
		     (= lep (line-end-position)) ; on initial line
		     ;; Not sure why it was like this -- fails in case of
		     ;; last internal function followed by first
		     ;; non-def statement of the main body.
                     ;;(and def-line (= in ci))
		     (= in ci)
		     (< in ci)))
	       (not (python-in-string/comment)))
	  (setq found t)))))

(defadvice python-beginning-of-defun
  (around python-beginning-of-defun-advice first (&rest rest) activate)
  ;; this strange incantation calls the original python function unless we're
  ;; in a derived mode.  The (setq ad-return-value ...) is how one modifies
  ;; the return value of advised functions.
  (if (not (pyrex-mode-p))
      ad-do-it
    (setq ad-return-value (apply 'pyrex-beginning-of-defun rest))))

(defun pyrex-end-of-defun ()
  "`end-of-defun-function' for Pyrex.
Finds end of innermost nested class or method definition."
  (let ((orig (point))
	(pattern (rx line-start (0+ space) (or "def" "cdef" "cpdef" "class") space)))
    ;; Go to start of current block and check whether it's at top
    ;; level.  If it is, and not a block start, look forward for
    ;; definition statement.
    (when (python-comment-line-p)
      (end-of-line)
      (forward-comment most-positive-fixnum))
    (if (not (python-open-block-statement-p))
	(python-beginning-of-block))
    (if (zerop (current-indentation))
	(unless (python-open-block-statement-p)
	  (while (and (re-search-forward pattern nil 'move)
		      (python-in-string/comment))) ; just loop
	  (unless (eobp)
	    (beginning-of-line)))
      ;; Don't move before top-level statement that would end defun.
      (end-of-line)
      (python-beginning-of-defun))
    ;; If we got to the start of buffer, look forward for
    ;; definition statement.
    (if (and (bobp) (not (looking-at (rx (or "def" "cdef" "cpdef" "class")))))
	(while (and (not (eobp))
		    (re-search-forward pattern nil 'move)
		    (python-in-string/comment)))) ; just loop
    ;; We're at a definition statement (or end-of-buffer).
    (unless (eobp)
      (python-end-of-block)
      ;; Count trailing space in defun (but not trailing comments).
      (skip-syntax-forward " >")
      (unless (eobp)			; e.g. missing final newline
	(beginning-of-line)))
    ;; Catch pathological cases like this, where the beginning-of-defun
    ;; skips to a definition we're not in:
    ;; if ...:
    ;;     ...
    ;; else:
    ;;     ...  # point here
    ;;     ...
    ;;     def ...
    (if (< (point) orig)
	(goto-char (point-max)))))

(defadvice python-end-of-defun
  (around python-end-of-defun-advice first (&rest rest) activate)
  ;; this strange incantation calls the original python function unless we're
  ;; in a derived mode.  The (setq ad-return-value ...) is how one modifies
  ;; the return value of advised functions.
  (if (not (pyrex-mode-p))
      ad-do-it
    (setq ad-return-value (apply 'pyrex-end-of-defun rest))))

;; Fixme: Consider top-level assignments, imports, &c.
(defun pyrex-current-defun ()
  "`add-log-current-defun-function' for Pyrex."
  (save-excursion
    ;; Move up the tree of nested `class' and `def' blocks until we
    ;; get to zero indentation, accumulating the defined names.
    (let ((start t)
	  accum)
      (while (or start (> (current-indentation) 0))
	(setq start nil)
	(python-beginning-of-block)
	(end-of-line)
	(beginning-of-defun)
	(if (looking-at (rx (0+ space) (or "def" "cdef" "cpdef" "class") (1+ space)
			    (group (1+ (or word (syntax symbol))))))
	    (push (match-string 1) accum)))
      (if accum (mapconcat 'identity accum ".")))))

(defadvice python-current-defun
  (around python-current-defun-advice first (&rest rest) activate)
  ;; this strange incantation calls the original python function unless we're
  ;; in a derived mode.  The (setq ad-return-value ...) is how one modifies
  ;; the return value of advised functions.
  (if (not (pyrex-mode-p))
      ad-do-it
    (setq ad-return-value (apply 'pyrex-current-defun rest))))

(provide 'pyrex)