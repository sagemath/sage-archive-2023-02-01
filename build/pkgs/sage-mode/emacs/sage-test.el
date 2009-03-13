;;; sage-test.el --- Test Sage files

;; Copyright (C) 2008  Nicholas Alexander

;; Author: Nicholas Alexander <ncalexan@pv109055.reshsg.uci.edu>
;; Keywords: sage test

(require 'sage)
(require 'sage-mode)

;; History of sage-test commands.
;;;###autoload
(defvar sage-test-history nil)

(defun sage-test-process-setup ()
  "Setup compilation variables and buffer for `sage-test'.
Set up `compilation-exit-message-function' and run `sage-test-setup-hook'."
  (set (make-local-variable 'compilation-exit-message-function)
       (lambda (status code msg)
	 (if (eq status 'exit)
	     (cond ((zerop code)
		    '("finished (all test passed)\n" . "all tests passed"))
		   ((> code 0)
		    '("finished with failing tests\n" . "failing tests"))
		   (t
		    (cons msg code)))
	   (cons msg code))))
  (run-hooks 'sage-test-setup-hook))

(defvar sage-test-regexp-alist
  '(("File \"\\(.*?\\)\", line \\([0-9]+\\):"
     1 2)
    ("File \"\\(.*?\\)\", line \\([0-9]+\\),"
     1 2 nil 1)
    ))

(defun sage-kill-compilation ()
  "Work hard to really kill sage-test."
  (interactive)
  (dotimes (dummy 50)
    ;; (get-buffer-process (current-buffer))
    (kill-compilation)
    (sleep-for 0 1)))

(define-compilation-mode sage-test-mode "sage-test"
  "Sets `grep-last-buffer' and `compilation-window-height'."
;;   (set (make-local-variable 'compilation-error-face)
;;        grep-hit-face)
;;   (set (make-local-variable 'compilation-disable-input) t)
  (local-set-key [(control c) (control k)] 'sage-kill-compilation)

  (set (make-local-variable 'compilation-error-regexp-alist)
       sage-test-regexp-alist)
  (set (make-local-variable 'compilation-process-setup-function)
       'sage-test-process-setup))

(defun sage-default-test-files ()
  (if (sage-mode-p)
      (buffer-file-name)
    (format "%s*" default-directory)))

(defun sage-default-test-command ()
  "Compute the default sage test command for `sage-test' to offer."
  (format "%s >/dev/null && %s -tp 4 %s" (sage-default-build-command) sage-command (sage-default-test-files)))

(defun sage-default-test-new-command ()
  "Compute the default sage test new command for `sage-test' to offer."
  (format "%s >/dev/null && %s -tnew" (sage-default-build-command) sage-command))

;;;###autoload
(defun sage-test (command-args)
  "Run sage-test, with user-specified args, and collect output in a buffer.
While sage-test runs asynchronously, you can use \\[next-error] (M-x next-error), or
\\<sage-test-mode-map>\\[compile-goto-error] in the sage-test
output buffer, to go to the lines where sage-test found matches.

With prefix arg, default to sage -tnew.

This command uses a special history list for its COMMAND-ARGS, so you can
easily repeat a sage-test command."
  (interactive
   (progn
     (let ((default (sage-default-test-command)))
       (list (read-from-minibuffer "Run sage-test (like this): "
				   (if current-prefix-arg (sage-default-test-new-command) default)
				   nil nil 'sage-test-history
				   (if current-prefix-arg nil default))))))

  ;; Setting process-setup-function makes exit-message-function work
  ;; even when async processes aren't supported.
  (compilation-start command-args 'sage-test-mode))

;;;_* "Interactive" doctesting
(defun sage-send-doctest-line-and-forward (&optional noshow)
  "If looking at a 'sage:' prompt, send this line and move to the next prompt
in this docstring.

If NOSHOW is nil, display the Sage process buffer."
  (interactive)
;;   (unless (python-in-string/comment)
;;     (error "Not in a Sage docstring"))
  (save-excursion
    (beginning-of-line)
    (unless (looking-at sage-test-prompt)
      (error "Not at a sage: prompt"))
    ;; send current line
    (re-search-forward sage-test-prompt)
    (sage-send-command (buffer-substring-no-properties
			(point) (line-end-position)) t)
    (unless noshow (display-buffer sage-buffer)))
  (save-restriction
    (narrow-to-defun)
    (end-of-line)
    (unless (re-search-forward sage-test-prompt (point-max) t)
      (forward-line 1))
    (end-of-line)))

(defun sage-send-all-doctest-lines (&optional noshow nogo)
  "If in a docstring, send every 'sage:' prompt.

If NOSHOW is nil, display the Sage process buffer.
If NOGO is nil, pop to the Sage process buffer."
  (interactive)
  (unless (python-in-string/comment)
    (error "Not in a Sage docstring"))
  (python-beginning-of-string)
  (re-search-forward sage-test-prompt)
  (ignore-errors
    (while t
      (sage-send-doctest-line-and-forward noshow)
      (inferior-sage-wait-for-prompt)))
  (unless nogo
    (pop-to-buffer sage-buffer)))

(defun sage-send-all-doctest-lines-in-file (&optional noshow nogo)
  "Go to the beginning of the file and send every 'sage:' prompt.

If NOSHOW is nil, display the Sage process buffer.
If NOGO is nil, pop to the Sage process buffer."
  (interactive)
  (save-restriction
    (goto-char (point-min))
    (ignore-errors
      (while t
	(re-search-forward sage-test-prompt)
	(sage-send-all-doctest-lines noshow t)
	(inferior-sage-wait-for-prompt))))
  (unless nogo
    (pop-to-buffer sage-buffer)))

(defun sage-retest-failing-files-from-buffer (buffer)
  (interactive "b")
  (save-match-data
    (with-current-buffer buffer
      (let ((files nil)
	    (root (concat (sage-root) "/")))
	(goto-char 0)
	(re-search-forward "^The following tests failed:")
	(while (re-search-forward "\\(devel/.*\\)" (point-max) t)
	  (setq files (cons (concat root (match-string 1)) files)))
	files))))

(defun sage-retest (buffer)
  "Retest only failing files scraped from an existing *sage-test* buffer.
"
  (interactive "bRetest failing files from buffer: ")
  (let* ((files (sage-retest-failing-files-from-buffer buffer))
	 (files-str (mapconcat #'identity files " ")))
    (flet ((sage-default-test-files () files-str))
      ;; (setq default-directory (sage-current-devel-root))
      (call-interactively 'sage-test))))

;;;###autoload
(defun sage-send-doctest (&optional all)
  "If looking at a sage: prompt, send the current doctest line to the inferior sage.
With prefix argument, send all doctests (at sage: prompts) until
the end of the docstring."
  (interactive "P")
  (if all
      (let ((current-prefix-arg nil))
	(sage-send-all-doctest-lines))
    (sage-send-doctest-line-and-forward)))

(provide 'sage-test)