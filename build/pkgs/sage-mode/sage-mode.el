;;; sage-mode.el --- Major mode for editing SAGE code

;; Copyright (C) 2007  Nick Alexander

;; Author: Nick Alexander <ncalexander@gmail.com>
;; Keywords: sage ipython python math

;; This file is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; This file is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to
;; the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
;; Boston, MA 02110-1301, USA.

;;; Commentary:

;;; Code:

(require 'python)
(require 'ansi-color)

;;; Use icicles for completing-read if possible
(require 'icicles nil t)

;;; Inferior SAGE major mode

(define-derived-mode
  inferior-sage-mode
  inferior-python-mode
  "Inferior SAGE"
  "Major mode for interacting with an inferior SAGE process."

  ; (setq sage-command (expand-file-name "~/bin/sage"))
  (setq comint-prompt-regexp
	(rx line-start (1+ (and (or "sage:" "....." ">>>" "...") " "))))
  (setq comint-redirect-finished-regexp "sage:") ; comint-prompt-regexp)
  ; ansi color doesn't play well with redirect
  ; (ansi-color-for-comint-mode-on)
  ; XXX what should be done here?
  ; (setq python-buffer sage-buffer)

  (define-key inferior-sage-mode-map
    [(control h) (control f)] 'ipython-describe-symbol)
)

(defcustom sage-command (expand-file-name "~/bin/sage")
  "Actual command used to run SAGE.
Additional arguments are added when the command is used by `run-sage' et al."
  :group 'sage
  :type 'string)

(defvar sage-buffer nil
  "*The current SAGE process buffer.

Commands that send text from source buffers to SAGE processes have
to choose a process to send to.  This is determined by buffer-local
value of `sage-buffer'.  If its value in the current buffer,
i.e. both any local value and the default one, is nil, `run-sage'
and commands that send to the Python process will start a new process.

Whenever \\[run-sage] starts a new process, it resets the default
value of `sage-buffer' to be the new process's buffer and sets the
buffer-local value similarly if the current buffer is in SAGE mode
or Inferior SAGE mode, so that source buffer stays associated with a
specific sub-process.

Use \\[sage-set-proc] to set the default value from a buffer with a
local value.")
(make-variable-buffer-local 'sage-buffer)

;;;###autoload
(defun run-sage (&optional cmd noshow new)
  "Run an inferior SAGE process, input and output via buffer *SAGE*.
CMD is the SAGE command to run.  NOSHOW non-nil means don't show the
buffer automatically.

Normally, if there is a process already running in `sage-buffer',
switch to that buffer.  Interactively, a prefix arg allows you to edit
the initial command line (default is `sage-command'); `-i' etc. args
will be added to this as appropriate.  A new process is started if:
one isn't running attached to `sage-buffer', or interactively the
default `sage-command', or argument NEW is non-nil.  See also the
documentation for `sage-buffer'.

Runs the hook `inferior-sage-mode-hook' \(after the
`comint-mode-hook' is run).  \(Type \\[describe-mode] in the process
buffer for a list of commands.)"
  (interactive (if current-prefix-arg
		   (list (read-string "Run SAGE: " sage-command) nil t)
		 (list sage-command)))
  (unless cmd (setq cmd sage-command))
  (setq sage-command cmd)
  ;; Fixme: Consider making `sage-buffer' buffer-local as a buffer
  ;; (not a name) in SAGE buffers from which `run-sage' &c is
  ;; invoked.  Would support multiple processes better.
  (let ((create-new-sage-p
	 (or new			; if you ask for it
	     (null sage-buffer)		; or there isn't a running sage
	     (not (comint-check-proc sage-buffer)) ; or there is a sage
					; buffer, but it's dead
	     )))
    (when create-new-sage-p
      (with-current-buffer
	  (let* ((cmdlist (python-args-to-list cmd))
		 ;; Set PYTHONPATH to import module emacs from emacs.py,
		 ;; but ensure that a user specified PYTHONPATH will
		 ;; override our setting, so that emacs.py can be
		 ;; customized easily.
		 (orig-path (getenv "PYTHONPATH"))
		 (path-sep (if (and orig-path (length orig-path)) ":" ""))
		 (data-path (concat "PYTHONPATH=" orig-path path-sep data-directory))
		 (process-environment
		  (cons data-path process-environment)))
	    (apply 'make-comint-in-buffer "SAGE"
		   (if new (generate-new-buffer "*SAGE*") "*SAGE*")
		   (car cmdlist) nil (cdr cmdlist)))
	;; Show progress
	(unless noshow (pop-to-buffer (current-buffer)))
	;; Update default SAGE buffers
	(setq-default sage-buffer (current-buffer))
	;; Update python-buffer too, so that evaluation keys work
	(setq-default python-buffer (current-buffer))
	;; Set up sensible prompt defaults, etc
	(inferior-sage-mode)
	(accept-process-output (get-buffer-process sage-buffer) 5)
	;; Ensure we're at a prompt before loading the functions we use
	;; XXX error-checking
	(python-send-receive-multiline "import emacs"))))

  ;; If we're coming from a sage-mode buffer, update inferior buffer
  (when (derived-mode-p 'sage-mode)
      (setq sage-buffer (default-value 'sage-buffer)) ; buffer-local
      ;; Update python-buffer too, so that evaluation keys work
      (setq python-buffer (default-value 'sage-buffer))) ; buffer-local

  ;; No matter how we got here, we want this inferior buffer to be the master
  ;; (when (comint-check-proc sage-buffer)
  ;;  (setq-default sage-buffer sage-buffer)
  ;;  (setq-default python-buffer sage-buffer))
  ;; Without this, help output goes into the inferior python buffer if
  ;; the process isn't already running.
  (sit-for 1 t)        ;Should we use accept-process-output instead?  --Stef
  (unless noshow (pop-to-buffer sage-buffer)))

;;; SAGE major mode

(provide 'sage)

(define-derived-mode
  sage-mode
  python-mode
  "SAGE"
  "Major mode for editing SAGE files."

  (define-key sage-mode-map [(control h) (control f)] 'ipython-describe-symbol)
)

;;;###autoload
(add-to-list 'interpreter-mode-alist '("sage" . sage-mode))
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.sage\\'" . sage-mode))

(defun python-qualified-module-name (file)
  "Find the qualified module name for filename FILE.

This recurses down the directory tree as long as there are __init__.py
files there, signalling that we are inside a package.

Returns a list of two elements.  The first is the top level package
directory; the second is the dotted Python module name.

Adapted from a patch posted to the python-mode trac."
  (let ((rec #'(lambda (d f)
		 (let* ((dir (file-name-directory d))
			(initpy (concat dir "__init__.py")))
		   (if (file-exists-p initpy)
		       (let ((d2 (directory-file-name d)))
			 (funcall rec (file-name-directory d2)
				  (concat (file-name-nondirectory d2) "." f)))
		     (list d f))))))
    (funcall rec (file-name-directory file)
	     (file-name-sans-extension (file-name-nondirectory file)))))

;;; Replace original `python-load-file' to use xreload and packages.
(ad-activate 'python-load-file)
(defadvice python-load-file
  (around nca-python-load-file first (file-name &optional xreload))
  "Load a Python file FILE-NAME into the inferior Python process.

THIS REPLACES THE ORIGINAL `python-load-file'.

If the file has extension `.py' import or reload it as a module.
Treating it as a module keeps the global namespace clean, provides
function location information for debugging, and supports users of
module-qualified names."
  (interactive
   (append (comint-get-source "Load Python file: " python-prev-dir/file
					   python-source-modes
					   t)
	   current-prefix-arg))	; because execfile needs exact name
  (comint-check-source file-name)     ; Check to see if buffer needs saving.
  (setq python-prev-dir/file (cons (file-name-directory file-name)
				   (file-name-nondirectory file-name)))
  (with-current-buffer (process-buffer (python-proc)) ;Runs python if needed.
    ;; Fixme: I'm not convinced by this logic from python-mode.el.
    (python-send-command
     (if (string-match "\\.py\\'" file-name)
	 (let* ((directory-module (python-qualified-module-name file-name))
		(directory (car directory-module))
		(module (cdr directory-module))
		(xreload-flag (if xreload "True" "False")))
	   (format "emacs.eimport(%S, %S, use_xreload=%s)"
		   module directory xreload-flag))
       (format "execfile(%S)" file-name)))
    (message "%s loaded" file-name)))

;;; Treat sage code as python source code
(add-to-list 'python-source-modes 'sage-mode)

(defun python-send-receive-multiline (command)
  "Send COMMAND to inferior Python (if any) and return result as a string.

This is an alternate `python-send-receive' that uses temporary buffers and
`comint-redirect-send-command-to-process'.
This implementation handles multi-line output strings gracefully.  At this
time, it does not handle multi-line input strings at all."
  (interactive "sCommand: ")
  (with-temp-buffer
    ;; Grab what Python has to say
    (comint-redirect-send-command-to-process
     command (current-buffer) (python-proc) nil t)
    ;; Wait for the redirection to complete
    (with-current-buffer (process-buffer (python-proc))
      (while (null comint-redirect-completed)
	(accept-process-output nil 1)))
    ;; Return the output
    (let ((output (buffer-substring-no-properties (point-min) (point-max))))
      (when (interactive-p)
	(message output))
      output)))

(defun python-current-word ()
  "Return python symbol at point."
  (with-syntax-table python-dotty-syntax-table
    (current-word)))

;;; `ipython-completing-read-symbol' is `completing-read' for python symbols
;;; using IPython's *? mechanism

(defvar ipython-completing-read-symbol-history ()
  "List of Python symbols recently queried.")

(defvar ipython-completing-read-symbol-pred nil
  "Default predicate for filtering queried Python symbols.")

(defvar ipython-completing-read-symbol-command "%s*?"
  "IPython command for generating completions.
Each completion should appear separated by whitespace.")

(defvar ipython-completing-read-symbol-cache ()
  "A list (last-queried-string string-completions).")

(defun ipython-completing-read-symbol-clear-cache ()
  "Clear the IPython completing read cache."
  (interactive)
  (setq ipython-completing-read-symbol-cache ()))

(defun ipython-completing-read-symbol-make-completions (string)
  "Query IPython for completions of STRING.

Return a list of completion strings.
Uses `ipython-completing-read-symbol-command' to query IPython."
  (let* ((command (format ipython-completing-read-symbol-command string))
	 (output (python-send-receive-multiline command)))
    (condition-case ()
	(split-string output)
      (error nil))))

(defun ipython-completing-read-symbol-function (string predicate action)
  "A `completing-read' programmable completion function for querying IPython.

See `try-completion' and `all-completions' for interface details."
  (let ((cached-string (first ipython-completing-read-symbol-cache))
	(completions (second ipython-completing-read-symbol-cache)))
    ;; Recompute table using IPython if neccessary
    (when (or (null completions)
	      (not (equal string cached-string)))
      (setq ipython-completing-read-symbol-cache
	    (list string (ipython-completing-read-symbol-make-completions string)))
      (setq completions (second ipython-completing-read-symbol-cache)))
    ;; Complete as necessary
    (if action
	(let ((all (all-completions string completions predicate)))
	  (if (eq action 'lambda)
	    (member string all)		; action is lambda
	    all))			; action is t
      (try-completion string completions predicate) ; action is nil
)))

(defun ipython-completing-read-symbol (&optional def require-match predicate)
  "Read a Python symbol (default: DEF) from user, completing with IPython.

Return a single element list, suitable for use in `interactive' forms.
If DEF is nil, default is `python-current-word'.
PREDICATE returns non-nil for potential completions.
See `completing-read' for REQUIRE-MATCH."
  (let* ((default (or def (python-current-word)))
	 (prompt (if (null default) "IPython symbol: "
		   (format "IPython symbol (default %s): " default)))
	 (func 'ipython-completing-read-symbol-function)
	 (pred (or predicate ipython-completing-read-symbol-pred))
	 (hist 'ipython-completing-read-symbol-history)
	 (enable-recursive-minibuffers t))
    (ipython-completing-read-symbol-clear-cache)
    (list (completing-read prompt func pred require-match nil hist default))))

;;; `ipython-describe-symbol' is `find-function' for python symbols using
;;; IPython's ? magic mechanism.

(defvar ipython-describe-symbol-not-found-regexp "Object `.*?` not found."
  "Regexp that matches IPython's 'symbol not found' warning.")

(defvar ipython-describe-symbol-command "%s?")

(defvar ipython-describe-symbol-temp-buffer-show-hook
  (lambda ()				; avoid xref stuff
    (toggle-read-only 1)
    (setq view-return-to-alist
	  (list (cons (selected-window) help-return-method))))
  "`temp-buffer-show-hook' installed for `ipython-describe-symbol' output.")

(defun ipython-describe-symbol-markup-function (string)
  "Markup IPython's inspection (?) for display."
  (when (string-match "[ \t\n]+\\'" string)
    (concat (substring string 0 (match-beginning 0)) "\n")))

(defun ipython-describe-symbol (symbol)
  "Get help on SYMBOL using IPython's inspection (?).
Interactively, prompt for SYMBOL."
  ;; Note that we do this in the inferior process, not a separate one, to
  ;; ensure the environment is appropriate.
  (interactive (ipython-completing-read-symbol nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (let* ((command (format ipython-describe-symbol-command symbol))
	 (raw-contents (python-send-receive-multiline command))
	 (help-contents
	  (or (ipython-describe-symbol-markup-function raw-contents)
	      raw-contents))
	 (temp-buffer-show-hook ipython-describe-symbol-temp-buffer-show-hook))
    ;; Handle symbol not found gracefully
    (when (string-match ipython-describe-symbol-not-found-regexp raw-contents)
      (error "Symbol not found"))
    (when (= 0 (length help-contents))
      (error "Symbol has no description"))
    ;; Ensure we have a suitable help buffer.
    (with-output-to-temp-buffer (help-buffer)
      (with-current-buffer standard-output
	;; Fixme: Is this actually useful?
	(help-setup-xref (list 'ipython-describe-symbol symbol) (interactive-p))
	(set (make-local-variable 'comint-redirect-subvert-readonly) t)
	(print-help-return-message)
	;; Finally, display help contents
	(princ help-contents)))))

;;; `sage-find-symbol' is `find-function' for SAGE.

(defun sage-find-symbol-command (symbol)
  "Return SAGE command to fetch position of SYMBOL."
  (format
   (concat "sage.misc.sageinspect.sage_getfile(%s), "
	   "sage.misc.sageinspect.sage_getsourcelines(%s)[-1]")
   symbol symbol))

(defvar sage-find-symbol-regexp "('\\(.*?\\)',[ \t\n]+\\([0-9]+\\))"
  "Match (FILENAME . LINE) from `sage-find-symbol-command'.")

(defun sage-find-symbol-noselect (symbol)
  "Return a pair (BUFFER . POINT) pointing to the definition of SYMBOL.

Queries SAGE to find the source file containing the definition of
FUNCTION in a buffer and the point of the definition.  The buffer
is not selected.

At this time, there is no error checking.  Later, if the function
definition can't be found in the buffer, returns (BUFFER)."
  (when (not symbol)
    (error "You didn't specify a symbol"))
  (let* ((command (sage-find-symbol-command symbol))
	 (raw-contents (python-send-receive-multiline command)))
    (unless (string-match sage-find-symbol-regexp raw-contents)
      (error "Symbol source not found"))
    (let ((filename (match-string 1 raw-contents))
	  (line-num (string-to-number (match-string 2 raw-contents))))
      (with-current-buffer (find-file-noselect filename)
	(goto-line line-num) ; XXX error checking?
	(cons (current-buffer) (point))))))

(defun sage-find-symbol-do-it (symbol switch-fn)
  "Find definition of SYMBOL in a buffer and display it.

SWITCH-FN is the function to call to display and select the
buffer."
      (let* ((orig-point (point))
	     (orig-buf (window-buffer))
	     (orig-buffers (buffer-list))
	     (buffer-point (save-excursion
			     (sage-find-symbol-noselect symbol)))
	     (new-buf (car buffer-point))
	     (new-point (cdr buffer-point)))
	(when buffer-point
	  (when (memq new-buf orig-buffers)
	    (push-mark orig-point))
	  (funcall switch-fn new-buf)
	  (when new-point (goto-char new-point))
	  (recenter find-function-recenter-line)
	  ;; (run-hooks 'find-function-after-hook)
	  )))

;;;###autoload
(defun sage-find-symbol (symbol)
  "Find the definition of the SYMBOL near point.

Finds the source file containing the defintion of the SYMBOL near point and
places point before the definition.
Set mark before moving, if the buffer already existed."
  (interactive (ipython-completing-read-symbol nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer))

;;;###autoload
(defun sage-find-symbol-other-window (symbol)
  "Find, in another window, the definition of SYMBOL near point.

See `sage-find-symbol' for details."
  (interactive (ipython-completing-read-symbol nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer-other-window))

;;;###autoload
(defun sage-find-symbol-other-frame (symbol)
  "Find, in another frame, the definition of SYMBOL near point.

See `sage-find-symbol' for details."
  (interactive (ipython-completing-read-symbol nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer-other-frame))
