;;; sage-build.el --- Build Sage

;; Copyright (C) 2008  Nicholas Alexander

;; Author: Nicholas Alexander <ncalexan@pv109055.reshsg.uci.edu>
;; Keywords: sage build

;;;###autoload
(defcustom sage-build-setup-hook nil
  "List of hook functions run by `sage-build-process-setup' (see `run-hooks')."
  :type 'hook
  :group 'sage-build)

;; History of sage-build commands.
;;;###autoload
(defvar sage-build-history nil)

(defun sage-build-process-setup ()
  "Setup compilation variables and buffer for `sage-build'.
Set up `compilation-exit-message-function' and run `sage-build-setup-hook'."
  (set (make-local-variable 'compilation-exit-message-function)
       (lambda (status code msg)
	 (if (eq status 'exit)
	     (cond ((zerop code)
		    '("finished (build failed)\n" . "build failed"))
		   ((> code 0)
		    '("finished (build succeeded)\n" . "build succeeded"))
		   (t
		    (cons msg code)))
	   (cons msg code))))
  (add-hook 'compilation-finish-functions
	    (lambda (buffer msg)
	      (when current-prefix-arg ;; XXX (and (string-match "succeed" msg) andrun)
		(run-sage t)
		)) nil t)
  (run-hooks 'sage-build-setup-hook))

(defvar sage-build-regexp-alist
  '(("File \"\\(.*?\\)\", line \\([0-9]+\\):"
     1 2)))

(define-compilation-mode sage-build-mode "sage-build"
  "Sets `grep-last-buffer' and `compilation-window-height'."
;;   (set (make-local-variable 'compilation-error-face)
;;        grep-hit-face)
;;   (set (make-local-variable 'compilation-disable-input) t)

  (set (make-local-variable 'compilation-error-regexp-alist)
       sage-build-regexp-alist)
  (set (make-local-variable 'compilation-process-setup-function)
       'sage-build-process-setup))

(defun sage-default-build-command ()
  "Compute the default sage build command for C-u M-x sage-build to offer."
  (format "%s -b" sage-command))

;;;###autoload
(defun sage-build (command-args)
  "Build sage (like sage -b), with user-specified args, and collect output in a buffer.
While sage-build runs asynchronously, you can use
\\[next-error] (M-x next-error), or
\\<sage-build-mode-map>\\[compile-goto-error] in the sage-build XXX
output buffer, to go to the lines where sage-build found matches.

With prefix arg (or andrun), act like sage -br: build sage and
spawn a new sage on success.

This command uses a special history list for its COMMAND-ARGS, so you can
easily repeat a sage-build command."
  (interactive
   (progn
     (let ((default (sage-default-build-command)))
       (list (read-from-minibuffer (if current-prefix-arg
				       "Build and run sage (like this): "
				     "Build sage (like this): ")
				    default
				    nil nil 'sage-build-history
				    (if current-prefix-arg nil default))))))

  ;; Setting process-setup-function makes exit-message-function work
  ;; even when async processes aren't supported.
  (compilation-start command-args 'sage-build-mode))

(provide 'sage-build)