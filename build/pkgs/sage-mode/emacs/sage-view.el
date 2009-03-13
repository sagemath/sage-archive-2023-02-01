;;; sage-view.el --- Typeset SAGE output on the fly
;; $Id: sage-view.el,v 1.2 2009/02/25 18:32:51 meulien Exp $

;; Copyright (C) 2008  Matthias Meulien, Nick Alexander

;; Authors: Matthias Meulien <matthias.meulien@xlim.fr>, Nick Alexander
;; <ncalexander@gmail.com>
;; Keywords: sage math image

;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; `sage-view' typesets output in an `inferior-sage-mode' buffer and displays
;; plots inline in an `inferior-sage-mode'.  Inline displays are context
;; sensitive; by default, right-clicking brings up a context-menu.

;; Use `sage-view' to enable the minor mode, and then
;; `sage-view-enable-inline-output', `sage-view-disable-inline-output' and
;; `sage-view-enable-inline-plots', `sage-view-disable-inline-plots' enable
;; and disable the relevant features.  You might add some of those functions
;; to `sage-view-hook' to configure `sage-view' to your liking.

;; You can customize `sage-view' using the Emacs customize interface by M-x
;; customize-group RET sage-view RET.  In particular you can customize
;; magnification and margins

;; This mode was inspired by doc-view.el by Tassilo Horn, preview.el
;; by David Kastrup, and imath.el by Yasuaki Honda.

;; The LaTex style used by preview.el is mandatory to use
;; sage-view.el. It is shipped with AUCTeX.

;;; Todo:
;; - Add a auto-reveal stuff to overlays
;; - Check that display is image capable
;; - Disabling sage-view mode should remove overlays
;; - Set color, center image, enlarge overlay to window full size
;; - Add zoom features to overlays
;; - Add horizontal scrolling

;; Bugs:
;; - Numpy output can be a text array... should not be inserted into
;;   $$ signs (hum... example?)

;;; Code:
(require 'sage)

(defvar sage-view-latex-head
  "\\documentclass{article}
\\usepackage[active, tightpage, pdftex, displaymath]{preview}
\\usepackage{amstext}
\\begin{document}
\\begin{preview}
\\begin{math}\n"
  "String to be inserted at the top of a LaTeX document.")

(defvar sage-view-latex-tail
  "\n\\end{math}
\\end{preview}
\\end{document}\n"
  "String to be inserted at the end of a LaTeX document.")

(defvar sage-view-start-string "<html><span class=\"math\">"
  "HTML tags that identify the begining of a math formula in Sage
  output.")

(defvar sage-view-final-string "</span></html>"
  "HTML tags that identify the end of a math formula in Sage
  output.")

(defvar sage-view-dir-name nil)
(defvar sage-view-inline-plots-enabled nil)
(defvar sage-view-inline-output-enabled nil)

(defun sage-view-latex->pdf (ov)
  "Start conversion of the LATEX document associated to OV to PDF.

See also `sage-view-process-overlay'."
  (let* ((latex (concat (overlay-get ov 'file-sans-extension) ".tex"))
	 (options (append
		   (list (concat "--output-directory=" sage-view-dir-name)
			 (concat "-interaction=" "nonstopmode")
			 (concat "-output-format=" "pdf")
			 latex)))
	 (proc (apply 'start-process
		      (append (list "latex->pdf" nil "latex") options))))
    (process-put proc 'overlay ov)
    (set-process-sentinel proc 'sage-view-latex->pdf-sentinel)))

(defun sage-view-pdf->png (ov)
  "Start conversion of the PDF file associated to OV to PNG.

See also `sage-view-process-overlay'."
  (let* ((base (overlay-get ov 'file-sans-extension))
	 (png (concat base ".png"))
	 (pdf (concat base ".pdf"))
	 (level (int-to-string sage-view-anti-aliasing-level))
	 (scale (or (overlay-get ov 'scale) sage-view-scale))
	 (options (append
		   sage-view-gs-options
		   (list (concat "-dTextAlphaBits=" level)
			 (concat "-dGraphicsAlphaBits=" level)
			 (concat "-sOutputFile=" png)
			 (concat "-r" (sage-view-compute-resolution scale))
			 pdf)))
	 (proc (apply 'start-process
		      (append
		       (list "pdf->png" "*sage-view*" sage-view-gs-command)
		       options))))
    (process-put proc 'overlay ov)
    (set-process-sentinel proc 'sage-view-pdf->png-sentinel)))

(defun sage-view-latex->pdf-sentinel (proc event)
  "If PROC (supposed to be a conversion process from LATEX to
PDF) was successful, convert the PDF to PNG.

See also `sage-view-process-overlay'."
  (let* ((ov (process-get proc 'overlay))
	 (base (overlay-get ov 'file-sans-extension)))
    (if (string-match "finished" event)
	(sage-view-pdf->png ov)
      (overlay-put ov 'display
		   (concat "Conversion failed (see " base ".log" ")")))))

(defun sage-view-pdf->png-sentinel (proc event)
  "If PROC (supposed to be a conversion process from PDF to PNG)
was successful, update the overlay associated to PROC.

See also `sage-view-process-overlay'."
  (let* ((ov (process-get proc 'overlay))
	 (base (overlay-get ov 'file-sans-extension))
	 (png (concat base ".png"))
	 (image (when (and (string-match "finished" event)
			   png (file-readable-p png))
		  (append (list 'image :type 'png :file png :margin sage-view-margin)))))
    (if (not image)
	(overlay-put ov 'display
		     (concat "Conversion failed (see " base ".log" ")"))
      (overlay-put ov 'display image))
    (sit-for 0)))

(defun sage-view-compute-resolution (scale)
  (let ((w (* scale (/ (* 25.4 (display-pixel-width))
	      (display-mm-width))))
	(h (* scale (/ (* 25.4 (display-pixel-height))
	      (display-mm-height)))))
    (concat (int-to-string w) "x" (int-to-string h))))

(defun sage-view-process-overlay (ov)
  "Associate a LATEX document to OV and start conversion process
from LATEX to PDF.

The conversion process is done by `sage-view-latex->pdf'. When it
ends, `sage-view-latex->pdf-sentinel' is called: If the
conversion is successful, a conversion process from PDF to PNG
starts. When it ends, `sage-view-pdf->png-sentinel' is called: If
the last conversion is successful, OV displays the resulting
image."
  (let* ((base (expand-file-name
		(make-temp-name "sage-view_") sage-view-dir-name))
	 (file (concat base ".tex")))
    (with-temp-file file
      (insert sage-view-latex-head)
      (insert (overlay-get ov 'math))
      (insert sage-view-latex-tail))
    (overlay-put ov 'file-sans-extension base)
    (sage-view-latex->pdf ov)))

(defun sage-view-output-filter-process-inline-output (string)
  "Substitute overlays to inline output.

Each region delimited by `sage-view-start-string' and
`sage-view-final-string' is replaced by an overlay.

This function expects the buffer to be narrowed to the current
output. And should be wrapped in a `save-excursion' and
`save-restriction' call.

See also `sage-view-output-filter'."
  (goto-char (point-min))
  (while (search-forward sage-view-start-string (point-max) t)
    (let* ((beg (point))
	   (end (- (search-forward sage-view-final-string (point-max) t)
		   (length sage-view-final-string)))
	   (text (buffer-substring-no-properties beg end))
	   (ov (make-overlay (- beg (length sage-view-start-string))
			     (+ end (length sage-view-final-string))
			     nil nil nil))
	   (map (make-sparse-keymap)))
      (overlay-put ov 'help-echo "mouse-3: Open contextual menu")
      (overlay-put ov 'math text)
      (define-key map [mouse-3]
	`(lambda (event) (interactive "e")
	   (sage-view-context-menu ,ov event)))
      (overlay-put ov 'keymap map)
      (sage-view-process-overlay ov))))

(defun sage-view-output-filter-process-inline-plots (string)
  "Generate and place one overlay image for one inline plot,
found by looking for a particular png file in directory
`sage-view-dir-name'.

This function expects the buffer to be narrowed to just the
current output; see `sage-view-output-filter' for how to do
that."
  ;; we need to foil emacs' image cache, which doesn't reload files with the same name
  (let* ((pngname (format "%s/sage-view.png" sage-view-dir-name))
	 (base (expand-file-name (make-temp-name "sage-view-plot_") sage-view-dir-name))
	 (pngname2 (concat base ".png")))
    ;; (message "Looking for plot at %s..." pngname)
    (if (not (and pngname
		  (file-exists-p pngname)
		  (file-readable-p pngname)))
	() ;; the not found branch
      ;; (message "Looking for plot at %s... not found." pngname)
      ;; the found branch
      ;; (message "Looking for plot at %s... found!" pngname)
      (dired-rename-file pngname pngname2 t)
      ;; (message "Renamed %s to %s" pngname pngname2)

      (goto-char (point-max))
      (let ((im (create-image pngname2 'png))) ;; adding a margin screws this up
	(setq sage-view-inline-plot-overlay
	      (make-overlay (- (point) 1) (- (point) 0) nil nil nil))
	(overlay-put sage-view-inline-plot-overlay 'display im)
	(overlay-put sage-view-inline-plot-overlay 'before-string "\n\n") ;; help alignment as much as possible
	(overlay-put sage-view-inline-plot-overlay 'after-string "\n\n")  ;; but emacs doesn't handle this right IMHO
      ;; (dired-delete-file pngname2 'always) ;; causes problems with emacs image loading
      ))))

(defun sage-view-output-filter (string)
  "Generate and place overlay images for inline output and inline plots.

Function to be inserted in `comint-output-filter-functions'."
  (save-excursion
    (save-restriction
      (narrow-to-region comint-last-input-end
			(process-mark (get-buffer-process (current-buffer))))
      (when sage-view-inline-plots-enabled
	(sage-view-output-filter-process-inline-plots string))
      (when (and sage-view-inline-output-enabled
		 (string-match (regexp-quote sage-view-final-string) string))
	(sage-view-output-filter-process-inline-output string)))))

;;;###autoload
(defun sage-view-enable-inline-output ()
  "Enable inline output pretty-printing, i.e. typeset output from sage in the `inferior-sage-mode' buffer.
WARNING: this communicates with the sage process.  Only use this
when `sage-view' mode is enabled and sage is running."
  (interactive)
  (setq sage-view-inline-output-enabled t)
  (python-send-receive-multiline "pretty_print_default(True);"))

;;;###autoload
(defun sage-view-disable-inline-output ()
  "Disable inline output pretty-printing, i.e. do not typeset output from sage in the `inferior-sage-mode' buffer.
WARNING: this communicates with the sage process.  Only use this
when `sage-view' mode is enabled and sage is running."
  (interactive)
  (setq sage-view-inline-output-enabled nil)
  (python-send-receive-multiline "pretty_print_default(False);"))

;;;###autoload
(defun sage-view-enable-inline-plots ()
  "Enable inline plotting, i.e. display plots in the `inferior-sage-mode' buffer and do not spawn an external viewer.
WARNING: this communicates with the sage process.  Only use this
when `sage-view' mode is enabled and sage is running."
  (interactive)
  (setq sage-view-inline-plots-enabled t)
  (python-send-receive-multiline "sage.plot.plot.DOCTEST_MODE = True;")
  (python-send-receive-multiline (format "sage.plot.plot.DOCTEST_MODE_FILE = '%s';"
					 (format "%s/sage-view.png" sage-view-dir-name))))

;;;###autoload
(defun sage-view-disable-inline-plots ()
  "Disable inline plotting, i.e. do not display plots in the `inferior-sage-mode' buffer and instead spawn an external viewer.
WARNING: this communicates with the sage process.  Only use this
when `sage-view' mode is enabled and sage is running."
  (interactive)
  (setq sage-view-inline-plots-enabled nil)
  (python-send-receive-multiline "sage.plot.plot.DOCTEST_MODE = False;")
  (python-send-receive-multiline (format "sage.plot.plot.DOCTEST_MODE_FILE = None;")))

(defun sage-view-create-temp ()
  "Create a temporary directory and set `sage-view-dir-name'
to its name.

Nothing is done if `sage-view-dir-name' equals the name of a
writable directory."
  (unless (and sage-view-dir-name
	       (file-directory-p sage-view-dir-name)
	       (file-writable-p sage-view-dir-name))
    (setq sage-view-dir-name
	  (make-temp-name (expand-file-name "tmp" "~/.sage/temp/")))
    (condition-case err
	(make-directory sage-view-dir-name)
      (error (message "Creation of `%s' failed: %s"
		      sage-view-dir-name (error-message-string err))))))

(defun sage-view-delete-temp ()
  "Delete the directory named after `sage-view-dir-name'."
  (when sage-view-dir-name
    (condition-case err
	(progn
	  (mapc #'delete-file
		(directory-files sage-view-dir-name t "^sage-view" t))
	  (delete-directory sage-view-dir-name))
      (error (message "Deletion of `%s' failed: %s"
		      sage-view-dir-name (error-message-string err))))))

(defun sage-view-overlay-activep (ov)
  "Check whether there is a valid image associated with OV."
  (eq (car (overlay-get ov 'display) 'image)))

(defun sage-view-copy-text (ov)
  "Copy LATEX source of OV into the kill buffer."
  (let ((text (overlay-get ov 'math)))
    (if text
	(kill-new text)
      (message "No LaTeX code available"))))

(defun sage-view-save-image (ov)
  "Copy image file associated to OV.

Make sure that there is a valid image associated with OV with
`sage-view-overlay-activep'."
  (let* ((spec (cdr (overlay-get ov 'display)))
	 (file (plist-get spec :file))
	 (name (when (and file (file-readable-p file))
		 (expand-file-name
		  (read-file-name "Write image to file: "
				 default-directory
				 "sage-view.png")))))
    (if name
	(copy-file file name))))

(defun sage-view-regenerate (ov)
  ""
  (overlay-put ov 'scale sage-view-scale)
  (sage-view-process-overlay ov))

(defun sage-view-zoom-in (ov)
  ""
  (let ((scale (or (overlay-get ov 'scale) sage-view-scale)))
    (overlay-put ov 'scale (+ scale sage-view-scale-factor))
    (message "Overlay's scale set to %s" scale)
    ;(overlay-put ov 'display (overlay-get ov 'math))
    (sage-view-process-overlay ov)))

(defun sage-view-zoom-out (ov)
  ""
  (let ((scale (or (overlay-get ov 'scale) sage-view-scale)))
    (overlay-put ov 'scale (- scale sage-view-scale-factor))
    (message "Overlay's scale set to %s" scale)
    ;(overlay-put ov 'display (overlay-get ov 'math))
    (sage-view-process-overlay ov)))

(defun sage-view-context-menu (ov ev)
  "Pop up a menu for OV at position EV."
  (popup-menu
   `("Sage Mode"
     ["Regenerate" (lambda () (interactive) (sage-view-regenerate ,ov))]
     ["Copy Text" (lambda () (interactive) (sage-view-copy-text ,ov))]
     ["Save As..." (lambda () (interactive) (sage-view-save-image ,ov))
      `(sage-view-overlay-activep ,ov)]
     ["Zoom in" (lambda () (interactive) (sage-view-zoom-in ,ov))
      `(sage-view-overlay-activep ,ov)]
     ["Zoom out" (lambda () (interactive) (sage-view-zoom-out ,ov))
      `(sage-view-overlay-activep ,ov)])
   ev))

;;;###autoload
(define-minor-mode sage-view
  "Toggle automatic typesetting of Sage output.

Typesetting of math formulas is done by LATEX subprocesses and
PDF to PNG conversions." nil
  :group 'sage-view
  :lighter " Sage-View"
  (if sage-view
      (progn
	(make-local-variable 'sage-view-dir-name)
	(sage-view-create-temp)
	(sage-view-enable-inline-output)
	(sage-view-enable-inline-plots)
	(make-local-variable 'comint-preoutput-filter-functions)
	(make-local-variable 'comint-output-filter-function)
	(add-hook 'comint-output-filter-functions 'sage-view-output-filter)
	(add-hook 'kill-buffer-hook 'sage-view-delete-temp))
    (progn
      (sage-view-disable-inline-output)
      (sage-view-disable-inline-plots)
      (remove-hook 'comint-output-filter-functions 'sage-view-output-filter)
      (remove-hook 'comint-preoutput-filter-functions 'sage-view-preoutput-filter))))

(provide 'sage-view)
;;; sage-view.el ends here
