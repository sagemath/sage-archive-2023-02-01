;;;_* sage-mode.el --- Major mode for editing SAGE code

;; Copyright (C) 2007, 2008  Nick Alexander

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

;; `sage-mode' is a major mode for editing sage (and python, and cython)
;; source code.  `inferior-sage-mode' is the companion mode for interacting
;; with a slave sage session.  See the help for `sage-mode' for help getting
;; started and the default key bindings.

;;; Code:

(require 'python)
(require 'comint)
(require 'ansi-color)
(require 'compile)

;;;_ + SAGE mode key bindings

(defvar sage-mode-map
  (let ((map (make-keymap)))	  ;`sparse' doesn't allow binding to charsets.
    (define-key map [(control c) (control c)] 'sage-send-buffer)
    (define-key map [(control c) (control r)] 'sage-send-region)
    (define-key map [(control c) (control j)] 'sage-send-doctest)
    (define-key map [(control c) (control t)] 'sage-test)
    (define-key map [(control c) (control b)] 'sage-build)
    (define-key map [(control h) (control f)] 'ipython-describe-symbol)
    (define-key map [(control h) (control g)] 'sage-find-symbol-other-window)

    (easy-menu-define menu-map map "Sage Mode menu"
      `("Sage"
	:help "sage-mode Specific Features"
	["Send Buffer" sage-send-buffer
	 :help "Send current buffer to inferior sage"]
	["Send Region" sage-send-region :active mark-active
	 :help "Send current region to inferior sage"]
	["Send Doctest" sage-send-doctest
	 :help "Send current doctest to inferior sage"]
	"-"
	["Run Sage" sage
	 :help "Run sage"]
	["Rerun Sage" sage-rerun
	 :help "Kill running sage and rerun sage"]
	["Build Sage" sage-build
	 :help "Build sage with \"sage -b\""]
	["Run Doctests" sage-test
	 :help "Run doctests with \"sage -t\""]))
    map)
  "Keymap for `sage-mode'.")

(defvar inferior-sage-mode-map
  (let ((map (make-keymap)))	  ;`sparse' doesn't allow binding to charsets.
    (define-key map [(control c) (control t)] 'sage-test)
    (define-key map [(control c) (control b)] 'sage-build)
    (define-key map [(control h) (control f)] 'ipython-describe-symbol)
    (define-key map [(control h) (control g)] 'sage-find-symbol-other-window)
    (define-key map (kbd "TAB") 'sage-pcomplete-or-help)

    (easy-menu-define menu-map map "Inferior Sage Mode menu"
      `("Sage"
	["Run Sage" sage
	 :help "Run sage"]
	["Rerun Sage" sage-rerun
	 :help "Kill running sage and rerun sage"]
	["Build Sage" sage-build
	 :help "Build sage with \"sage -b\""]
	["Run Doctests" sage-test
	 :help "Run doctests with \"sage -t\""]))
    map)
  "Keymap for `inferior-sage-mode'.")

;;;_* Inferior SAGE major mode for interacting with a slave SAGE process

;;;###autoload
(define-derived-mode
  inferior-sage-mode
  inferior-python-mode
  "Inferior SAGE"
  "Major mode for interacting with an inferior SAGE process."

  (sage-set-buffer (current-buffer))

  (setq comint-prompt-regexp inferior-sage-prompt)
  (setq comint-prompt-read-only t)
  (setq comint-redirect-finished-regexp comint-prompt-regexp)

  (setq comint-input-sender 'sage-input-sender)
  ;; I type x? a lot
  (set 'comint-input-filter 'sage-input-filter)

  (make-local-variable 'compilation-error-regexp-alist)
  (make-local-variable 'compilation-error-regexp-alist-alist)
  (add-to-list 'compilation-error-regexp-alist-alist sage-test-compilation-regexp)
  (add-to-list 'compilation-error-regexp-alist 'sage-test-compilation)
  (add-to-list 'compilation-error-regexp-alist-alist sage-build-compilation-regexp)
  (add-to-list 'compilation-error-regexp-alist 'sage-build-compilation)

  (pcomplete-sage-setup)
  (compilation-shell-minor-mode 1))

(defun inferior-sage-wait-for-prompt ()
  "Wait until the SAGE process is ready for input."
  (with-current-buffer sage-buffer
    (let* ((sprocess (get-buffer-process sage-buffer))
	    (success nil)
	    (timeout 0))
      (while (progn
	       (if (not (eq (process-status sprocess) 'run))
		   (error "SAGE process has died unexpectedly.")
		 (if (> (setq timeout (1+ timeout)) inferior-sage-timeout)
		     (error "Timeout waiting for SAGE prompt. Check inferior-sage-timeout."))
		 (accept-process-output nil 0 1)
		 (sit-for 0)
		 (goto-char (point-max))
		 (forward-line 0)
		 (setq success (looking-at inferior-sage-prompt))
		 (not (or success (looking-at ".*\\?\\s *"))))))
      (goto-char (point-max))
      success)))

(defun sage-input-filter (string)
  "A `comint-input-filter' that keeps all input in the history."
  t)

;;;_* IPython magic commands

(defcustom ipython-input-handle-magic-p t
  "Non-nil means handle IPython magic commands specially."
  :group 'ipython
  :type 'boolean)

(defvar ipython-input-string-is-magic-regexp
  "\\(\\**\\?\\??\\)\\'"
  "Regexp matching IPython magic input.

The first match group is used to dispatch handlers in
`ipython-input-handle-magic'.")

(defun ipython-input-string-is-magic-p (string)
  "Return non-nil if STRING is IPython magic."
  (string-match ipython-input-string-is-magic-regexp string))

(defvar ipython-input-magic-handlers '(("**?" . ipython-handle-magic-**?)
				       ("??"  . ipython-handle-magic-??)
				       ("?"  . ipython-handle-magic-?))
  "Association list (STRING . FUNCTION) of IPython magic handlers.

Each FUNCTION should take arguments (PROC STRING MATCH) and
return non-nil if magic input was handled, nil if input should be
sent normally.")

;; 			  (cons nil (cdr apropos-item)))))
;; 	  (insert-text-button (symbol-name symbol)
;; 			      'type 'apropos-symbol
;; 			      ;; Can't use default, since user may have
;; 			      ;; changed the variable!
;; 			      ;; Just say `no' to variables containing faces!
;; 			      'face apropos-symbol-face)

(define-button-type 'sage-apropos-command
  'face 'bold
  'apropos-label "Command:"
  'help-echo "mouse-2, RET: Display more help on this command"
  'follow-link t
  'action (lambda (button)
	    (ipython-describe-symbol (button-get button 'apropos-symbol))))

(defun sage-apropos-mode ()
  (interactive)
  (apropos-mode)
  (save-excursion
    (goto-char (point-min))
    (let ((case-fold-search nil))
      (while (re-search-forward "^ \\* \\(.*?\\) \\(.*?\\):" nil t)
	(let ((type (match-string 1))
	      (name (match-string 2)))
	  ;; reformat everything
	  (replace-match "")
	  (insert
	   (format "%s\n  " (propertize name 'face 'bold)))
	  (insert-text-button (format "%s:" type)
			      'type 'sage-apropos-command
			      'apropos-label (format "%s:" type)
			      'face apropos-label-face
			      'apropos-symbol name)
	  )))))

(defun sage-apropos (symbol)
  (interactive "sApropos SAGE command: ")
  (when (or (null symbol) (equal "" symbol))
    (error "No command"))
  (with-output-to-temp-buffer "*sage-apropos*"
   (with-current-buffer "*sage-apropos*"
     (python-send-receive-to-buffer (format "apropos('%s')" symbol) "*sage-apropos*")
     (sage-apropos-mode)
     (goto-char 0)
     (insert (format "Sage apropos for %s:\n\n" symbol))
     t)))

(defun ipython-handle-magic-**? (proc string &optional match)
  "Handle SAGE apropos **?."
  (when (string-match "\\(.*?\\)\\*\\*\\?" string)
    (sage-apropos (match-string 1 string))))

(defun ipython-handle-magic-? (proc string &optional match)
  "Handle IPython magic ?."
  (when (string-match "\\(.*?\\)\\?" string)
    (ipython-describe-symbol (match-string 1 string))))

(defun ipython-handle-magic-?? (proc string &optional match)
  "Handle IPython magic ??."
  (when (string-match "\\(.*?\\)\\?\\?" string)
    (sage-find-symbol-other-window (match-string 1 string))))

(defun ipython-input-handle-magic (proc string)
  "Handle IPython magic input STRING in process PROC.

Return non-nil if input was handled; nil if input should be sent
normally."
  (when (string-match ipython-input-string-is-magic-regexp string)
    (let* ((match (match-string 1 string))
	   (handler (cdr (assoc match ipython-input-magic-handlers))))
      (when handler
	(condition-case ()
	    ;; I can't explain why, but this seems to work perfectly with ??
	    (save-excursion
	      (funcall handler proc string match))
	  ;; XXX print error message?
	  (error nil))))))

(defun ipython-input-sender (proc string)
  "Function for sending to process PROC input STRING.

When `ipython-input-handle-magic-p' is non-nil, this uses
`ipython-input-string-is-magic-p' to look for ipython magic
commands, such as %prun, etc, and magic suffixes, such as ? and
??, and handles them... magically?  It hands them off to
`ipython-input-handle-magic' for special treatment.

Otherwise, `comint-simple-send' just sends STRING plus a newline."
  (if (and ipython-input-handle-magic-p ; must want to handle magic
	       (ipython-input-string-is-magic-p string) ; ... be magic
	       (ipython-input-handle-magic proc string)) ; and be handled
      ;; To have just an input history creating, clearing new line entered
      (comint-simple-send proc "")
    (comint-simple-send proc string)))	; otherwise, you're sent

;;;_* Make it easy to read long output
(defun sage-input-sender (proc string)
  "A `comint-input-sender' that might send output to a buffer."
  (when (not current-prefix-arg)
    ;; When not prefixed, send normally
    (ipython-input-sender proc string))
  (when current-prefix-arg
    ;; When prefixed, output to buffer -- modelled after shell-command
    (let ((output-buffer (get-buffer-create "*Sage Output*")))
      (save-excursion
	;; Clear old output -- maybe this is a bad idea?
	(set-buffer output-buffer)
	(setq buffer-read-only nil)
	(erase-buffer))
      (python-send-receive-to-buffer string output-buffer)
      (comint-simple-send proc "") ;; Clears input, saves history
      (if (with-current-buffer output-buffer (> (point-max) (point-min)))
	  (display-message-or-buffer output-buffer)))))

;;;_* SAGE process management

(defun sage-send-startup-command ()
  (sage-send-command sage-startup-command t))
(add-hook 'sage-startup-hook 'sage-send-startup-command)

(defvaralias 'sage-buffer 'python-buffer)
;; (defvar sage-buffer nil
;;   "*The current SAGE process buffer.

;; Commands that send text from source buffers to SAGE processes have
;; to choose a process to send to.  This is determined by buffer-local
;; value of `sage-buffer'.  If its value in the current buffer,
;; i.e. both any local value and the default one, is nil, `run-sage'
;; and commands that send to the Python process will start a new process.

;; Whenever \\[run-sage] starts a new process, it resets the default
;; value of `sage-buffer' to be the new process's buffer and sets the
;; buffer-local value similarly if the current buffer is in SAGE mode
;; or Inferior SAGE mode, so that source buffer stays associated with a
;; specific sub-process.

;; Use \\[sage-set-proc] to set the default value from a buffer with a
;; local value.")
;; (make-variable-buffer-local 'sage-buffer)

(defun sage-mode-p ()
  (interactive)
  (derived-mode-p 'sage-mode 'pyrex-mode))

(defun inferior-sage-mode-p ()
  (interactive)
  (derived-mode-p 'inferior-sage-mode))

(defun sage-all-inferior-sage-buffers ()
  "List of names of sage buffers."
  (let ((bufs nil))
    (save-excursion
      (dolist (buf (buffer-list))
	(with-current-buffer buf
	  (when (inferior-sage-mode-p)
	    (push (buffer-name buf) bufs))))
      (nreverse bufs))))

(defun sage-running-inferior-sage-buffers ()
  (let ((bufs nil))
    (save-excursion
      (dolist (buf (sage-all-inferior-sage-buffers))
	(with-current-buffer buf
	  (when (get-buffer-process (current-buffer))
	    (push (current-buffer) bufs))))
      bufs)))

(defun sage-mode-line-name-for-sage-buffer (buffer)
  (format ": [%s]" (buffer-name buffer)))

(defun sage-update-mode-line (buffer)
  (setq mode-line-process (sage-mode-line-name-for-sage-buffer buffer))
  (force-mode-line-update))

(defalias 'set-sage-buffer 'sage-set-buffer)
(defun sage-set-buffer (buffer)
  (interactive
   (list (progn
	   ;; (unless (sage-mode-p)
;; 	     (error "Not in a sage-mode buffer!"))
	   (completing-read
	    "SAGE buffer: " (sage-all-inferior-sage-buffers) nil nil
	    (car (sage-all-inferior-sage-buffers))))))
  (let ((chosen-buffer (with-current-buffer buffer (current-buffer))))
    (setq sage-buffer chosen-buffer)
    (setq python-buffer chosen-buffer) ; update python-buffer too
    (when (sage-mode-p)
      (sage-update-mode-line chosen-buffer))))

(defun python-proc ()
  "Return the current Python process.
See variable `python-buffer'.  Starts a new process if necessary."
  ;; Fixme: Maybe should look for another active process if there
  ;; isn't one for `python-buffer'.

  (cond ((inferior-sage-mode-p)
	 ;; if we're in a sage buffer, that's us
	 (get-buffer-process (current-buffer)))
	((comint-check-proc sage-buffer)
	 ;; if we refer to a good sage instance, use it
	 (get-buffer-process sage-buffer))
	((sage-running-inferior-sage-buffers)
	 ;; if there are other good sage instances, use one of them
	 (call-interactively 'sage-set-buffer)
	 (python-proc))
	(t
	 ;; otherwise, start a new sage and try again
	 (run-sage nil t)
	 (python-proc))))

;; History of sage-run commands.
;;;###autoload
(defvar sage-run-history nil)

(defun sage-create-new-sage (cmd)
  (interactive
   (progn
     (let ((default (or cmd sage-command)))
       (list (read-from-minibuffer "Run sage (like this): "
				   default
				   nil nil 'sage-run-history
				   default)))))
  (unless cmd
    (setq cmd sage-command))
  (with-current-buffer
      (let* ((sage-buffer-base-name (format "*SAGE-%s*" (sage-current-branch)))
	     (sage-buffer-name (if new (generate-new-buffer sage-buffer-base-name) sage-buffer-base-name))
	     (cmdlist (python-args-to-list cmd))
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
	       sage-buffer-name
	       (car cmdlist) nil (cdr cmdlist)))))

(defun sage-new-sage-p ()
  (interactive)
  (or (and (inferior-sage-mode-p)	; in a sage buffer but it's dead
	   (not (comint-check-proc (current-buffer))))
      (or (null sage-buffer)		      ; if there isn't a running sage
	  (not (comint-check-proc sage-buffer))))) ; or the sage buffer is dead

;;;###autoload
(defalias 'sage 'run-sage)
;;;###autoload
(defalias 'sage-run 'run-sage)
;;;###autoload
(defun run-sage (&optional new cmd noshow)
  "Run an inferior SAGE process, input and output via buffer *SAGE*.

NEW non-nil means always create a new buffer and SAGE process.
CMD is the SAGE command to run.
NOSHOW non-nil means don't show the buffer automatically.

Normally, if there is a process already running in `sage-buffer',
switch to that buffer.  A new process is started if: one isn't
running attached to `sage-buffer', or interactively the default
`sage-command', or argument NEW is non-nil.  See also the
documentation for `sage-buffer'.

Runs the hook `inferior-sage-mode-hook' \(after the
`comint-mode-hook' is run).  \(Type \\[describe-mode] in the process
buffer for a list of commands.)"
  (interactive "P")
  ;; Fixme: Consider making `sage-buffer' buffer-local as a buffer
  ;; (not a name) in SAGE buffers from which `run-sage' &c is
  ;; invoked.  Would support multiple processes better.
  (if (not (or new (sage-new-sage-p)))
      (unless noshow (pop-to-buffer sage-buffer))
    (setq sage-buffer (if (called-interactively-p)
			  (call-interactively 'sage-create-new-sage)
			(sage-create-new-sage cmd)))
    (set-default 'sage-buffer sage-buffer) ; update defauls
    (set-default 'python-buffer python-buffer)

    (with-current-buffer sage-buffer
      (unless noshow (pop-to-buffer sage-buffer)) ; show progress
      (unless (inferior-sage-mode-p)
	(inferior-sage-mode))
      (when (inferior-sage-wait-for-prompt) ; wait for prompt
	(run-hooks 'sage-startup-hook)))

    (when (sage-mode-p)
      ;; If we're coming from a sage-mode buffer, update inferior buffer
      (message "Buffer %s will use sage %s" (current-buffer) sage-buffer)
      (sage-set-buffer sage-buffer))
    (unless noshow (pop-to-buffer sage-buffer))))

(defun sage-set-buffer-name ()
  (interactive)
  "Change the current SAGE buffer name to include the current branch."
  (when (sage-current-branch)
    (rename-buffer
     (generate-new-buffer-name (format "*SAGE-%s*" (sage-current-branch))))))

(defun sage-root ()
  "Return SAGE_ROOT."
  (interactive)
  (save-match-data
    (let ((lst (split-string (shell-command-to-string (concat sage-command " -root")))))
    (nth 0 lst))))

(defun sage-current-branch-link ()
  "Return the current SAGE branch link, i.e., the target of devel/sage."
  (interactive)
  (save-match-data
    (let ((lst (split-string (shell-command-to-string (concat sage-command " -branch")))))
      (if (= 1 (length lst))
	  (nth 0 lst)
	"main"))))

(defun sage-current-branch ()
  "Return the current SAGE branch name."
  (interactive)
  (save-match-data
    (if (and (inferior-sage-mode-p)
	     (string-match "\\*SAGE-\\(.*\\)\\*" (buffer-name)))
	(match-string 1 (buffer-name))
      (sage-current-branch-link))))

(defun sage-current-devel-root ()
  (interactive)
  "Return the current SAGE branch directory."
  (format "%s/devel/sage-%s" (sage-root) (sage-current-branch)))

;;;_* SAGE major mode for editing SAGE library code

;;;###autoload
(define-derived-mode
  sage-mode
  python-mode
  "SAGE"
  "Major mode for editing SAGE files.

The major entry points are:

`sage', to spawn a new sage session.

`sage-send-buffer', to send the current buffer to the inferior sage, using
\"load\"; `sage-send-region', to send the current region to the inferior
sage, using \"load\"; and `sage-send-doctest', to send the docstring point is
currently looking at to the inferior sage interactively.

`sage-test', to execute \"sage -t\" and friends and parse the output

`sage-build', to execute \"sage -b\" and friends and parse the output.

`sage-rerun' to restart an inferior sage in an existing buffer, and
`sage-build' with a prefix argument to execute \"sage -br\" to rebuild sage
and restart a fresh inferior sage in an existing buffer.

\\{sage-mode-map}"
  (set (make-local-variable 'font-lock-multiline) t)
  (set (make-local-variable 'font-lock-defaults)
       `(, (cons
	    (cons "XXX\\(.*\n\\)*?XXX" font-lock-comment-face)
	    python-font-lock-keywords)
	 nil nil nil nil
	 (font-lock-syntactic-keywords . python-font-lock-syntactic-keywords)
	 ))
)

(defun sage-font-lock ()
  "Install Sage font-lock patterns."
  (interactive)
  ;; (font-lock-add-keywords 'sage-mode python-font-lock-keywords 'set) ;; XXX
;;   (font-lock-add-keywords 'sage-mode
;; 			  `(("\\(\\*\\*\\)test\\(\\*\\*\\)" . 'font-lock-comment-face)))
)

(add-hook 'sage-mode-hook 'sage-font-lock)

;;;_* Treat SAGE code as Python source code

;;;###autoload
(add-to-list 'interpreter-mode-alist '("sage" . sage-mode))
;;;###autoload
(add-to-list 'auto-mode-alist '("\\.sage\\'" . sage-mode))

(add-to-list 'python-source-modes 'sage-mode)

;;;###autoload
(defun sage-send-buffer ()
  "Send the current buffer to the inferior sage process.
The buffer is loaded using sage's \"load\" command."
  (interactive)
  (when (buffer-file-name)
    ;; named file -- offer to save it, then send it
    (when (buffer-modified-p)
      (save-some-buffers))
    (sage-send-command (format "load %s" (buffer-file-name)) t))
  (unless (buffer-file-name)
    ;; un-named buffer -- use sage-send-region
    (sage-send-region (point-min) (point-max)))
  (pop-to-buffer sage-buffer))

;;;###autoload
(defun sage-send-region (start end)
  "Send the region to the inferior sage process.
The region is treated as a temporary \".sage\" file with minimal
processing.  The logic is that this command is intended to
emulate interactive input, although this isn't perfect: sending
the region \"2\" does not print \"2\"."
  ;; The region is evaluated from a temporary file.  This avoids
  ;; problems with blank lines, which have different semantics
  ;; interactively and in files.  It also saves the inferior process
  ;; buffer filling up with interpreter prompts.  We need a Python
  ;; function to remove the temporary file when it has been evaluated
  ;; (though we could probably do it in Lisp with a Comint output
  ;; filter).  This function also catches exceptions and truncates
  ;; tracebacks not to mention the frame of the function itself.
  ;;
  ;; The `compilation-shell-minor-mode' parsing takes care of relating
  ;; the reference to the temporary file to the source.
  ;;
  ;; Fixme: Write a `coding' header to the temp file if the region is
  ;; non-ASCII.
  (interactive "r")
  (let* ((f (make-temp-file "sage" nil ".sage"))
	 (command (format "load '%s' # loading region..." f))
	 (orig-start (copy-marker start)))
    (when (save-excursion
	    (goto-char start)
	    (/= 0 (current-indentation))) ; need dummy block
      (save-excursion
	(goto-char orig-start)
	;; Wrong if we had indented code at buffer start.
	(set-marker orig-start (line-beginning-position 0)))
      (write-region "if True:\n" nil f nil 'nomsg))
    (write-region start end f t 'nomsg)
    (message "Sending region to sage buffer...")
    (sage-send-command command t) ;; the true is whether to show the input line or not
    (pop-to-buffer sage-buffer)
    (message "Sending region to sage buffer... done.")
    (with-current-buffer (process-buffer (python-proc))
      ;; Tell compile.el to redirect error locations in file `f' to
      ;; positions past marker `orig-start'.  It has to be done *after*
      ;; `python-send-command''s call to `compilation-forget-errors'.
      (compilation-fake-loc orig-start f))))

;;;_* Integrate SAGE mode with Emacs

;;;###autoload
(defun sage-pcomplete-or-help ()
  "If point is after ?, describe preceding symbol; otherwise, pcomplete."
  (interactive)
  (if (not (looking-back "[^\\?]\\?"))
      (pcomplete)
    (save-excursion
      (backward-char)
      (when (python-current-word)
	(ipython-describe-symbol (python-current-word))))))

;;;_ + Set better grep defaults for SAGE and Pyrex code

;; (eval-after-load "grep"
;;   (progn
;;     (add-to-list 'grep-files-aliases '("py" . "{*.py,*.pyx}"))
;;     (add-to-list 'grep-files-aliases '("pyx" . "{*.py,*.pyx}"))))

;;;_ + Make devel/sage files play nicely, and don't jump into site-packages if possible

;;; It's annoying to get lost in sage/.../site-packages version of files when
;;; `sage-find-symbol' and friends jump to files.  It's even more annoying when
;;; the file is not correctly recognized as sage source!

(add-to-list 'auto-mode-alist '("devel/sage.*?\\.py\\'" . sage-mode))
(add-to-list 'auto-mode-alist '("devel/sage.*?\\.pyx\\'" . pyrex-mode))

(defvar sage-site-packages-regexp "\\(local.*?site-packages.*?\\)/sage"
  "Regexp to match sage site-packages files.

Match group 1 will be replaced with devel/sage-branch")

(add-hook 'find-file-hook 'sage-warn-if-site-packages-file)
(defun sage-warn-if-site-packages-file()
  "Warn if sage FILE is in site-packages and offer to find current branch version."
  (let ((f (buffer-file-name (current-buffer))))
    (and f (string-match sage-site-packages-regexp f)
         (if (y-or-n-p "This is a sage site-packages file, open the real file? ")
             (sage-jump-to-development-version)
           (push '(:propertize "SAGE-SITE-PACKAGES-FILE:" face font-lock-warning-face)
                 mode-line-buffer-identification)))))

(defun sage-development-version (filename)
  "If FILENAME is in site-packages, current branch version, else FILENAME."
  (save-match-data
    (let* ((match (string-match sage-site-packages-regexp filename)))
      (if (and filename match)
	  ;; handle current branch somewhat intelligiently
	  (let* ((base (concat (substring filename 0 (match-beginning 1)) "devel/"))
		 (branch (or (file-symlink-p (concat base "sage")) "sage")))
	    (concat base branch (substring filename (match-end 1))))
	filename))))

(defun sage-jump-to-development-version ()
  "Jump to current branch version of current FILE if we're in site-packages version."
  (interactive)
  (let ((filename (sage-development-version (buffer-file-name (current-buffer))))
	(maybe-buf (find-buffer-visiting filename)))
    (if maybe-buf (pop-to-buffer maybe-buf)
      (find-alternate-file filename))))

(require 'advice)
(defadvice compilation-find-file
  (before sage-compilation-find-file (marker filename directory &rest formats))
  "Always try to find compilation errors in FILENAME in the current branch version."
  (ad-set-arg 1 (sage-development-version filename)))
(ad-activate 'compilation-find-file)

;;;_ + Integrate eshell with hg

(defadvice hg-root
  (before eshell-hg-root (&optional path))
  "Use current directory in eshell-mode for hg-root if possible.
Use current devel directory in inferior-sage-mode for hg-root if possible."
  (when (derived-mode-p 'eshell-mode)	; buffer local in eshell buffers
    (ad-set-arg 0 default-directory))
  (when (inferior-sage-mode-p) ; buffer local in inferior sage buffers
    (ad-set-arg 0 (sage-current-devel-root))))

(ad-activate 'hg-root)

;;;_ + Integrate with eshell

(defconst sage-test-compilation-regexp
  (list 'sage-test-compilation
	"^File \"\\(.*\\)\", line \\([0-9]+\\)"
	1
	2))

(defconst sage-build-compilation-regexp
  (list 'sage-build-compilation
	"^\\(.*\\):\\([0-9]+\\):\\([0-9]+\\):"
	1 2 3))

;; To add support for SAGE build and test errors to *compilation* buffers by
;; default, evaluate the following four lines.
;;
;; (add-to-list 'compilation-error-regexp-alist-alist sage-test-compilation-regexp)
;; (add-to-list 'compilation-error-regexp-alist 'sage-test-compilation)
;; (add-to-list 'compilation-error-regexp-alist-alist sage-build-compilation-regexp)
;; (add-to-list 'compilation-error-regexp-alist 'sage-build-compilation)

(defun eshell-sage-command-hook (command args)
  "Handle some SAGE invocations specially.

Without ARGS, run SAGE in an emacs `sage-mode' buffer.

With first ARGS starting with \"-b\" or \"-t\", run SAGE in an
emacs `compilation-mode' buffer.

Otherwise (for example, with ARGS \"-hg\", run SAGE at the eshell
prompt as normal.

This is an `eshell-named-command-hook' because only some parameters modify the
command; other times, it has to execute as a standard eshell command."
  (when (equal command "sage")
    (cond ((not args)
	   ;; run sage inside emacs
	   (run-sage nil sage-command nil)
	   t)
	  ((member (substring (car args) 0 2) '("-t" "-b"))
	   ;; echo sage build into compilation buffer
	   (throw 'eshell-replace-command
		  (eshell-parse-command
		   "compile"
		   (cons "sage" (eshell-flatten-list args))))))))
(add-hook 'eshell-named-command-hook 'eshell-sage-command-hook)

;; From http://www.emacswiki.org/cgi-bin/wiki?EshellFunctions
(defun eshell/compile (&rest args)
  "Use `compile' to do background makes."
  (if (eshell-interactive-output-p)
      (let ((compilation-process-setup-function
	     (list 'lambda nil
		   (list 'setq 'process-environment
			 (list 'quote (eshell-copy-environment))))))
	(compile (eshell-flatten-and-stringify args))
	(pop-to-buffer compilation-last-buffer))
    (throw 'eshell-replace-command
	   (let ((l (eshell-stringify-list (eshell-flatten-list args))))
	     (eshell-parse-command (car l) (cdr l))))))
(put 'eshell/compile 'eshell-no-numeric-conversions t)

;;;_* Load relative modules correctly

(defun python-qualified-module-name (file)
  "Find the qualified module name for filename FILE.

This recurses down the directory tree as long as there are __init__.py
files there, signalling that we are inside a package.

Returns a pair (PACKAGE . MODULE).  The first is the top level
package directory; the second is the dotted Python module name.

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
(defadvice python-load-file
  (around nca-python-load-file first (file-name &optional no-xreload))
  "Load a Python file FILE-NAME into the inferior Python process.

Without prefix argument, use fancy new xreload. With prefix
argument, use default Python reload.

THIS REPLACES THE ORIGINAL `python-load-file'.

If the file has extension `.py' import or reload it as a module.
Treating it as a module keeps the global namespace clean, provides
function location information for debugging, and supports users of
module-qualified names."
  (interactive
   (append (comint-get-source
	    (format "%s Python file: " (if current-prefix-arg "reload" "xreload"))
	    python-prev-dir/file
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
		(xreload-flag (if no-xreload "False" "True")))
	   (format "emacs.eimport(%S, %S, use_xreload=%s)"
		   module directory xreload-flag))
       (format "execfile(%S)" file-name)))
    (message "%s loaded" file-name)))
(ad-activate 'python-load-file)

;;;_* Convenient *programmatic* Python interaction

(defvar python-default-tag-noerror "_XXX1XXX_NOERROR")
(defvar python-default-tag-error "_XXX1XXX_ERROR")

(defun python-protect-command (command &optional tag-noerror tag-error)
  "Wrap Python COMMAND in a try-except block and signal error conditions.

Print TAG-NOERROR on successful Python execution and TAG-ERROR on
error conditions."
  (let* ((tag-noerror (or tag-noerror python-default-tag-noerror))
	 (tag-error   (or tag-error   python-default-tag-error))
	 (lines (split-string command "\n"))
	 (indented-lines
	  (mapconcat (lambda (x) (concat "    " x)) lines "\n")))
    (format "try:
%s
    print %S,
except e:
    print e
    print %S,
" indented-lines tag-noerror tag-error)))

(defmacro with-python-output-to-buffer (buffer command &rest body)
  "Send COMMAND to inferior Python, redirect output to BUFFER, and execute
BODY in that buffer.

The value returned is the value of the last form in body.

Block while waiting for output."
  (declare (indent 2) (debug t))
  `(with-current-buffer ,buffer
     ;; Grab what Python has to say
     (comint-redirect-send-command-to-process
      (python-protect-command ,command)
      (current-buffer) (python-proc) nil t)
     ;; Wait for the redirection to complete
     (with-current-buffer (process-buffer (python-proc))
       (while (null comint-redirect-completed)
	 (accept-process-output nil 1)))
     (message (buffer-name))
     ;; Execute BODY
     ,@body
     ))

(defmacro with-python-output-buffer (command &rest body)
  "Send COMMAND to inferior Python and execute BODY in temp buffer with
  output.

The value returned is the value of the last form in body.

Block while waiting for output."
  (declare (indent 1) (debug t))
  `(with-temp-buffer
     (with-python-output-to-buffer (current-buffer) ,command
       ,@body)))

;; (with-python-output-to-buffer "*scratch*" "x?\x?"
;;   (message "%s" (buffer-name)))

(defun sage-send-command (command &optional echo-input)
  "Evaluate COMMAND in inferior Python process.

If ECHO-INPUT is non-nil, echo input in process buffer."
  (interactive "sCommand: ")
  (with-current-buffer (process-buffer (python-proc))
    (goto-char (point-max))
    (if echo-input
	(with-current-buffer (process-buffer (python-proc))
	  ;; Insert and evaluate input string in place
	  (let ((old (comint-get-old-input-default)))
	    (delete-field)
	    (insert command)
	    (comint-send-input nil t)
	    (insert old)))
      (python-send-command command))))

(defun python-send-receive-to-buffer (command buffer &optional echo-output)
  "Send COMMAND to inferior Python (if any) and send output to BUFFER.

If ECHO-OUTPUT is non-nil, echo output to process buffer.

This is an alternate `python-send-receive' that uses temporary buffers and
`comint-redirect-send-command-to-process'.  Block while waiting for output.
This implementation handles multi-line output strings gracefully.  At this
time, it does not handle multi-line input strings at all."
  (interactive "sCommand: ")
  (with-current-buffer buffer
    ;; Grab what Python has to say
    (comint-redirect-send-command-to-process
     command (current-buffer) (python-proc) echo-output t)
    ;; Wait for the redirection to complete
    (with-current-buffer (process-buffer (python-proc))
      (while (null comint-redirect-completed)
	(accept-process-output nil 1)))))

(defun python-send-receive-multiline (command)
  "Send COMMAND to inferior Python (if any) and return result as a string.

This is an alternate `python-send-receive' that uses temporary buffers and
`comint-redirect-send-command-to-process'.  Block while waiting for output.
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

;;;_* Generally useful tidbits

(defun python-current-word ()
  "Return python symbol at point."
  (with-syntax-table python-dotty-syntax-table
    ;; the t makes current-word strict: returns nil if point is not in the word
    (current-word t)))

;;;_* IPython and SAGE completing reads

;;;_ + `ipython-completing-read-symbol' is `completing-read' for python symbols
;;; using IPython's *? mechanism

(defvar ipython-completing-read-symbol-history ()
  "List of Python symbols recently queried.")

(defvar ipython-completing-read-symbol-pred nil
  "Default predicate for filtering queried Python symbols.")

(defvar ipython-completing-read-symbol-command "_ip.IP.magic_psearch('-a %s*')"
  "IPython command for generating completions.
Each completion should appear separated by whitespace.")

(defvar ipython-completing-read-symbol-cache ()
  "A pair (LAST-QUERIED-STRING . COMPLETIONS).")

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
  (let ((cached-string (car ipython-completing-read-symbol-cache))
	(completions (cdr ipython-completing-read-symbol-cache)))
    ;; Recompute table using IPython if neccessary
    (when (or (null completions)
	      (not (equal string cached-string)))
      (setq ipython-completing-read-symbol-cache
	    (cons string (ipython-completing-read-symbol-make-completions string)))
      (setq completions
	    (cdr ipython-completing-read-symbol-cache)))
    ;; Complete as necessary
    (cond ((eq action 'lambda) ; action is 'lambda
	   (test-completion string completions))
	  (action   ; action is t
	   (pcomplete-uniqify-list (all-completions string completions predicate)))
	  (t	    ; action is nil
	   (try-completion string completions predicate)))))

(defun ipython-completing-read-symbol
  (&optional prompt def require-match predicate)
  "Read a Python symbol (default: DEF) from user, completing with IPython.

Return a single element list, suitable for use in `interactive' forms.
PROMPT is the prompt to display, without colon or space.
If DEF is nil, default is `python-current-word'.
PREDICATE returns non-nil for potential completions.
See `completing-read' for REQUIRE-MATCH."
  (let* ((default (or def (python-current-word)))
	 (prompt (if (null default) (concat prompt ": ")
		   (concat prompt (format " (default %s): " default))))
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

(define-button-type 'help-sage-function-def
  :supertype 'help-xref
  'help-function #'sage-find-symbol-other-window
  'help-echo (purecopy "mouse-2, RET: find function's definition"))

(defun ipython-describe-symbol-markup-buffer (symbol)
  "Markup IPython's inspection (?) in current buffer for display."
  (help-make-xrefs (current-buffer))
  (save-excursion
    (save-match-data
      (let ((case-fold-search nil))
	;; Make HEADERS: stand out
	(goto-char (point-min))
	(while (re-search-forward "\\([A-Z][^a-z]+\\):" nil t) ;; t means no error
	  (toggle-read-only 0)
	  (add-text-properties (match-beginning 1) (match-end 1) '(face bold)))

	;; make File: a link
	(goto-char (point-min))
	(while (re-search-forward "File:\\s-*\\(.*\\)" nil t) ;; t means no error
	  (toggle-read-only 0)
	  (replace-match (sage-development-version (match-string 1)) nil nil nil 1)
	  (help-xref-button 1 'help-sage-function-def symbol)
	  (toggle-read-only 1))
	t))))

(defun ipython-describe-symbol (symbol)
  "Get help on SYMBOL using IPython's inspection (?).
Interactively, prompt for SYMBOL."
  ;; Note that we do this in the inferior process, not a separate one, to
  ;; ensure the environment is appropriate.
  (interactive (ipython-completing-read-symbol "Describe symbol" nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (let* ((command (format ipython-describe-symbol-command symbol))
	 (raw-contents (python-send-receive-multiline command))
	 (help-contents
	  (or (ipython-describe-symbol-markup-function raw-contents)
	      raw-contents))
	 (temp-buffer-show-hook ipython-describe-symbol-temp-buffer-show-hook))
    ;; XXX Handle exceptions; perhaps (with-python-output ...) or similar
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
	(princ help-contents)))
    ;; Markup help buffer
    (with-current-buffer (help-buffer)
      (ipython-describe-symbol-markup-buffer symbol))))

;;;_ + `sage-find-symbol' is `find-function' for SAGE.

(defun sage-find-symbol-command (symbol)
  "Return SAGE command to fetch position of SYMBOL."
  (format
   (concat "sage.misc.sageinspect.sage_getfile(%s), "
	   "sage.misc.sageinspect.sage_getsourcelines(%s)[-1] + 1")
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
    (let* ((raw-filename (match-string 1 raw-contents))
	   (filename (sage-development-version raw-filename))
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
	  t)))

;;;###autoload
(defun sage-find-symbol (symbol)
  "Find the definition of the SYMBOL near point.

Finds the source file containing the defintion of the SYMBOL near point and
places point before the definition.
Set mark before moving, if the buffer already existed."
  (interactive (ipython-completing-read-symbol "Find symbol" nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer))

;;;###autoload
(defun sage-find-symbol-other-window (symbol)
  "Find, in another window, the definition of SYMBOL near point.

See `sage-find-symbol' for details."
  (interactive (ipython-completing-read-symbol "Find symbol" nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer-other-window))

;;;###autoload
(defun sage-find-symbol-other-frame (symbol)
  "Find, in another frame, the definition of SYMBOL near point.

See `sage-find-symbol' for details."
  (interactive (ipython-completing-read-symbol "Find symbol" nil t))
  (when (or (null symbol) (equal "" symbol))
    (error "No symbol"))
  (sage-find-symbol-do-it symbol 'switch-to-buffer-other-frame))

;;;_ + `try-complete-sage-symbol-partially' is a `hippie-expand' function for SAGE

(defun he-sage-symbol-beg ()
  (save-excursion
    (with-syntax-table python-dotty-syntax-table
      (skip-syntax-backward "w_")
      (point))))

(defun he-sage-symbol-end ()
  (save-excursion
    (with-syntax-table python-dotty-syntax-table
      (skip-syntax-forward "w_")
      (point))))

(defun try-complete-sage-symbol-partially (old)
  "Try to complete as a SAGE symbol, as many characters as unique.

The argument OLD is nil if this is the first call to this
function, non-nil if this is a subsequent call.

Returns t if a unique, possibly partial, completion is found; nil
otherwise."
  (let ((expansion nil))
    (when (not old)
      (he-init-string (he-sage-symbol-beg) (point))
      (unless (string= he-search-string "")
	(setq expansion (ipython-completing-read-symbol-function
			 he-search-string nil nil)))
      (when (or (eq expansion nil)
		(string= expansion he-search-string)
		(he-string-member expansion he-tried-table))
	(setq expansion nil)))
    (if (not expansion)
	(progn
	  (when old (he-reset-string))
	  nil)
        (progn
	  (he-substitute-string expansion)
	  t))))

;;;_ + `pcomplete' support

;; if pcomplete is available, set it up!
;; (when (featurep 'pcomplete)

(defun pcomplete-sage-setup ()
  (interactive)
  (set (make-local-variable 'pcomplete-autolist)
      nil)
  (set (make-local-variable 'pcomplete-cycle-completions)
      nil)
  (set (make-local-variable 'pcomplete-use-paring)
      nil)

  (set (make-variable-buffer-local 'pcomplete-default-completion-function)
       'pcomplete-sage-default-completion)
  (set (make-variable-buffer-local 'pcomplete-command-completion-function)
       'pcomplete-sage-default-completion)
  (set (make-variable-buffer-local 'pcomplete-parse-arguments-function)
       'pcomplete-parse-sage-arguments)

  (set (make-variable-buffer-local 'pcomplete-termination-string)
       "")
  )

(defun pcomplete-sage-completions ()
  (save-excursion
    (save-restriction
      (let ((stub (nth pcomplete-index pcomplete-args)))
	(when (and stub (not (string= stub "")))
	  (ipython-completing-read-symbol-clear-cache)
	  (ipython-completing-read-symbol-function stub nil t))))))

(defun pcomplete-sage-default-completion ()
  (pcomplete-here (pcomplete-sage-completions)))

(defun pcomplete-parse-sage-arguments ()
  (list
   (list (buffer-substring-no-properties (he-sage-symbol-beg)
					 (he-sage-symbol-end)))
   (he-sage-symbol-beg)))

;;;_* Make it easy to sagetest files and methods

(defun sage-test-file-inline (file-name &optional method)
  "Run sage-test on file FILE-NAME, with output to underlying the SAGE buffer.

We take pains to test the correct module.

If METHOD is non-nil, try to test only the single method named METHOD.
Interactively, try to find current method at point."
  (interactive
   (append
    (comint-get-source "Load SAGE file: "
		       python-prev-dir/file python-source-modes t))
   (list current-prefix-arg))
  (comint-check-source file-name)     ; Check to see if buffer needs saving.
  (setq python-prev-dir/file (cons (file-name-directory file-name)
				   (file-name-nondirectory file-name)))
  (let* ((directory-module (python-qualified-module-name file-name))
	 (directory (car directory-module))
	 (module (cdr directory-module))
	 (command (format "sage.misc.sagetest.sagetest(%s)" module)))
    (sage-send-command command nil)))

(defun sage-test-file-to-buffer (file-name &optional method)
  "Run sage-test on file FILE-NAME, with output to a new buffer.

We take pains to test the correct module.

If METHOD is non-nil, try to test only the single method named METHOD.
Interactively, try to find current method at point."
  (interactive
   (append
    (comint-get-source "Load SAGE file: "
		       python-prev-dir/file python-source-modes t))
   (list current-prefix-arg))
  (comint-check-source file-name)     ; Check to see if buffer needs saving.
  (setq python-prev-dir/file (cons (file-name-directory file-name)
				   (file-name-nondirectory file-name)))
  (let* ((directory-module (python-qualified-module-name file-name))
	 (directory (car directory-module))
	 (module (cdr directory-module))
	 (command (format "sage.misc.sagetest.sagetest(%s)" module))
	 (compilation-error-regexp-alist '(sage-test sage-build)))
    (with-temp-buffer
      (compile (eshell-flatten-and-stringify args))
      (python-send-receive-to-buffer command (current-buffer)))))

(defvar sage-test-file 'sage-test-file-to-buffer)

;;;_* Read Mercurial's .hg bundle files naturally

(define-derived-mode
  mercurial-bundle-mode
  fundamental-mode
  "Mercurial .hg bundle"
  "Major mode for interacting with an inferior SAGE process."
  (completing-read "Against repository: " '("sage-main" "sage-blah" "sage-nuts") nil t "sage-main")
  nil
)
(add-to-list 'auto-mode-alist '("\\.hg\\'" . mercurial-bundle-mode))

;;;_* Setup imenu by default
(when (featurep 'imenu)
  (add-hook 'sage-mode-hook 'imenu-add-menubar-index))

(provide 'sage-mode)
