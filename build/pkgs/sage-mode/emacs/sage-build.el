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

(defcustom sage-rerun-command (format "%s" sage-command)
  "Actual command used to rerun SAGE.
Additional arguments are added when the command is used by `rerun-sage' et al."
  :group 'sage
  :type 'string)

(defun sage-wait-until-dead (process seconds msg)
  (let ((killed nil))
    (message "Trying %s..." msg)
    (accept-process-output)
    (if (with-timeout (seconds (not killed))
	  (while (not killed)
	    (accept-process-output sprocess 1)
	    (setq killed (not (memq (process-status process) '(open run stop))))
	    ;; (setq killed (eq (process-status process) 'exit))

	    ;; (goto-char (point-max))
	    ;; (beginning-of-line)
	    ;; (setq killed (looking-at "Process SAGE killed"))
	    ))
	(progn
	  (accept-process-output sprocess 0 1)
	  (message "Trying %s... failed!" msg)
	  nil)
      (accept-process-output sprocess 0 1)
      (message "Trying %s... done." msg)
      t
      )))

(defalias 'restart-sage 'rerun-sage)

(defun rerun-sage ()
  (interactive)
  (let ((bufs (sage-all-inferior-sage-buffers))
	(buffer nil))
    (cond ((not bufs)
	   (progn
	     (message "Need a SAGE buffer to rerun SAGE in...")
	     (call-interactively 'run-sage)))
	  ((= 1 (length bufs))
	   (setq buffer (car bufs)))
	  (t (setq buffer
		   (completing-read
		    "Rerun SAGE in buffer: " bufs nil nil
		    (car bufs)))))
    ;; (kill-buffer buffer)
    (with-current-buffer buffer
      ;; (message (get-buffer-process buffer))
      (let* ((sprocess (get-buffer-process (current-buffer))))
	(when (and sprocess
		   (not (eq (process-status sprocess) 'exit)))
	  ;; (set-process-sentinel sprocess nil)
	  (comint-kill-input)
	  (comint-send-eof)
	  (sage-wait-until-dead sprocess 2 "soft kill")

	  (when (not (eq (process-status sprocess) 'exit))
	    (comint-kill-subjob)
	    (sage-wait-until-dead sprocess 3 "hard kill"))))

      ;; (comint-mode)
      (goto-char (point-max))
      (insert "\nRestarting SAGE...\n\n")
      (goto-char (point-max))
      (run-sage nil sage-rerun-command) ;; restart
      (goto-char (point-max))
      )))

	  ;; 	  ;; if we're not dead yet...
;; 	  (message "Killing soft")
;; 	  (goto-char (point-max))
;; 	  (comint-kill-input)
;; 	  (comint-send-eof) ;; kill nicely
;; 	  (accept-process-output nil 0 1)
;; 	  (sit-for 0.5) ;; XXX correct way to do this?
;; 	  (if (not (eq (process-status sprocess) 'exit))
;; 	      (message "Killing hard")
;; 	      (comint-kill-subjob) ;; kill rudely
;; 	    (accept-process-output nil 0 1))
;; 	  (sit-for 0.5)) ;; XXX correct way to do this?

;;       (pop-to-buffer buffer)
;; ;;       (while (and (get-buffer-process (current-buffer))
;; ;; 		  (< dummy 10))
;; ;; 	(comint-send-eof) ;; nicer than comint-kill-subjob
;; ;; 	(sleep-for 0 1)
;; ;; 	(incf dummy))
;;       (message "A")
;;       (when (get-buffer-process buffer)
;; 	(with-current-buffer buffer
;; 	  (comint-kill-subjob)))
;;       (run-sage t)
;;       (unless (get-buffer-process buffer)
;; 	(kill-buffer buffer))
;;       (message "C"))))

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
		(rerun-sage)
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