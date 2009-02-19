;;; sage-view.el --- Typeset SAGE output on the fly

;; Copyright (C) 2008  Matthias Meulien

;; Author: Matthias Meulien <matthias.meulien@xlim.fr>
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

;; Put this file in a directory in `load-path' and add the following
;; line to `.emacs':

;; if you want sage-view to be enabled for every sage session, try
;; (add-hook 'sage-startup-hook 'sage-view-always)

;; This mode was inspired by doc-view.el by Tassilo Horn, preview.el
;; by David Kastrup, and imath.el by Yasuaki Honda.

;; The LaTex style used by preview.el is mandatory to use
;; sage-view.el. It is shipped with AUCTeX.

;;; Todo:
;; - Add a auto-reveal stuff to overlays so that one can copy-cut from
;;   them (use convertion from pdf to text?)
;; - Check that display is image capable
;; - Disabling sage preview should remove overlays
;; - Set color, center image, enlarge overlay to window full size
;; - Split output string according to <html> tags (split-string)
;; - Add zoom features to overlays
;; - Add XML parser error treatment
;; - Add horizontal scrolling
;; - Make variables local, so that one can run multiple instance of SAGE
;; - Delete files

;; Bugs:
;; - Seems that latex produce buggy DVI when there are overlays... It
;;   happens with large matrices of polynomials
;; - error in process filter: Wrong type argument: number-or-marker-p, nil
;; - Multiple output
;; - Numpy output can be a text array... should not be inserted into $$
;;   signs

;;; Code:
(require 'sage)

;; (add-hook 'inferior-sage-mode-hook 'sage-view)

(defvar sage-view-head
  "\\documentclass{article}\\usepackage[active, tightpage, pdftex, displaymath]{preview}\\usepackage{amstext}\\begin{document}\\begin{preview}\$")
;; it should be possible to choose the conversion technique, and
;; change pdftex to dvips
(defvar sage-view-tail
  "\$\\end{preview}\\end{document}\n")

(defvar sage-view-temp-dir-name nil
  "Name of directory for temporary files")

(defvar sage-view-conversion-process nil)
(defvar sage-view-conversion-buffer "*sage-view*")

(defvar sage-view-anti-aliasing-level 2)
(defvar sage-view-ghostscript-program "gs")
(defvar sage-view-ghostscript-options
  (list "-sDEVICE=png16m"
	(concat "-dTextAlphaBits=" (int-to-string sage-view-anti-aliasing-level))
	"-dBATCH"
	(concat "-dGraphicsAlphaBits=" (int-to-string sage-view-anti-aliasing-level))
	"-dSAFER"
	"-q"
	"-dNOPAUSE"))

;; XXX these should all be customizable
(defvar sage-view-resolution nil)
(defvar sage-view-scale 1.2)
(defvar sage-view-current-overlay nil)

(defvar sage-view-start-string "<html><span class=\"math\">")
(defvar sage-view-final-string "</span></html>")

(defvar sage-view-inline-plots-enabled nil) ;; do we want inline plotting?
(defvar sage-view-inline-output-enabled nil) ;; do we want inline latex output?

(defun sage-view-latex->dvi (latex)
  "Convert LATEX to DVI asynchronously."
  (setq sage-view-conversion-process
	(apply 'start-process
	       (append (list "latex->dvi" sage-view-conversion-buffer
			     "latex")
		       (list (concat "--output-directory=" (shell-quote-argument (sage-view-temp-dir)))
			     (concat "-interaction=" (shell-quote-argument "nonstopmode"))
			     (concat "-output-format=" (shell-quote-argument "dvi")))
		       (list latex))))
  (set-process-sentinel sage-view-conversion-process
			'sage-view-latex->dvi-sentinel)
  (process-put sage-view-conversion-process 'dvi-file
	       (concat (file-name-sans-extension latex) ".dvi")))

(defun sage-view-latex->pdf (latex)
  "Convert LATEX to PDF asynchronously."
  (setq sage-view-conversion-process
	(apply 'start-process
	       (append (list "latex->pdf" sage-view-conversion-buffer
			     "latex")
		       (list (concat "--output-directory=" (shell-quote-argument (sage-view-temp-dir)))
			     (concat "-interaction=" (shell-quote-argument "nonstopmode"))
			     (concat "-output-format=" (shell-quote-argument "pdf")))
		       (list latex))))
  (set-process-sentinel sage-view-conversion-process
			'sage-view-latex->pdf-sentinel)
  (process-put sage-view-conversion-process 'pdf-file
	       (concat (file-name-sans-extension latex) ".pdf")))

(defun sage-view-dvi->ps (dvi ps)
  "Convert DVI to PS asynchronously."
  (setq sage-view-conversion-process
	(start-process "dvi->ps" sage-view-conversion-buffer
		       "dvips" "-Pwww" "-A" "-o" ps dvi))
  (set-process-sentinel sage-view-conversion-process
			'sage-view-dvi->ps-sentinel)
  (process-put sage-view-conversion-process 'ps-file ps))


(defun sage-view-pdf/ps->png (ps-pdf png)
  "Convert PDF-PS to PNG asynchronously."
  (setq sage-view-conversion-process
	(apply 'start-process
	       (append (list "pdf/ps->png" sage-view-conversion-buffer
			     sage-view-ghostscript-program)
		       sage-view-ghostscript-options
		       (list (concat "-sOutputFile=" png))
		       (list (concat "-r" (sage-view-compute-resolution)))
		       (list ps-pdf))))
  (set-process-sentinel sage-view-conversion-process
			'sage-view-pdf/ps->png-sentinel)
  (process-put sage-view-conversion-process 'png-file png))

(defun sage-view-latex->dvi-sentinel (proc event)
  "If LATEX->DVI conversion was successful, convert the DVI to PS."
  (let* ((dvi (process-get proc 'dvi-file))
	 (ps (concat (file-name-sans-extension dvi) ".dvi")))
    (if (string-match "finished" event)
	(sage-view-dvi->ps dvi ps)
      (overlay-put sage-view-current-overlay 'display
		   (concat "SAGE View failed (sage-view-latex->dvi-sentinel) (see "
			   (file-name-sans-extension dvi) ".log" ")")))))

(defun sage-view-latex->pdf-sentinel (proc event)
  "If LATEX->PDF conversion was successful, convert the PDF to PNG."
  (let* ((pdf (process-get proc 'pdf-file))
	 (png (concat (file-name-sans-extension pdf) ".png")))
    (if (string-match "finished" event)
	(sage-view-pdf/ps->png pdf png)
      (overlay-put sage-view-current-overlay 'display
		   (concat "SAGE View failed (sage-view-latex->pdf-sentinel) (see "
			   (file-name-sans-extension pdf) ".log" ")")))))

(defun sage-view-dvi->ps-sentinel (proc event)
  "If DVI->PS conversion was successful, convert the PS to PNG."
  (let* ((ps (process-get proc 'ps-file))
	 (png (concat (file-name-sans-extension ps) ".png")))
    (if (string-match "finished" event)
	(sage-view-pdf/ps->png ps png)
      (overlay-put sage-view-current-overlay 'display
		   (concat "SAGE View failed (sage-view-dvi->ps-sentinel) (see "
			   (file-name-sans-extension ps) ".log" ")")))))

(defun sage-view-pdf/ps->png-sentinel (proc event)
  "If PDF/PS->PNG conversion was successful, update
  `sage-view-current-overlay' overlay."
  (let ((png (process-get proc 'png-file)))
    (if (string-match "finished" event)
	(let ((image (if (and png (file-readable-p png))
			(create-image png 'png))))
	  (cond
	   (image
	    (overlay-put sage-view-current-overlay 'display image)
	    (sit-for 0))
	   (t (overlay-put sage-view-current-overlay 'display
			   (concat "SAGE View failed (sage-view-pdf/ps->png-sentinel A) (see "
			 (file-name-sans-extension png) ".log" ")")))))
      (overlay-put sage-view-current-overlay 'display
		   (concat "SAGE View failed (sage-view-pdf/ps->png-sentinel B) (see "
			 (file-name-sans-extension png) ".log" ")")))))

(defun sage-view-compute-resolution ()
  (let ((w (* sage-view-scale (/ (* 25.4 (display-pixel-width))
	      (display-mm-width))))
	(h (* sage-view-scale (/ (* 25.4 (display-pixel-height))
	      (display-mm-height)))))
    (concat (int-to-string w) "x" (int-to-string h))))

(defun sage-view-output-filter-process-inline-output (string)
  "Generate and place overlay images for inline output delimited
by `sage-view-start-string' and `sage-view-final-string'.

This function expects the buffer to be narrowed to just the
current output; see `sage-view-output-filter' for how to do
that."

  (goto-char (point-min))

  (when (string-match (regexp-quote sage-view-final-string) string) ;; we found the end of some html; go back and texify that range
    (setq sage-view-text (buffer-substring-no-properties (point-min) (point-max)))

    (when (search-forward sage-view-start-string (point-max) t) ;; change this to while to do many
      (setq sage-view-overlay-start (- (point) (length sage-view-start-string)))
      (search-forward sage-view-final-string) ;; find the terminal
      (search-backward sage-view-final-string) ;; position ourselves at the front of the terminal
      (setq sage-view-overlay-final (+ (point) (length sage-view-final-string)))
      (setq sage-view-insert
	    (buffer-substring-no-properties
	     (+ sage-view-overlay-start (length sage-view-start-string))
	     (- sage-view-overlay-final (length sage-view-final-string))))
      ;; (message "sage-view-insert _%s_" sage-view-insert)

      (setq sage-view-current-overlay (make-overlay sage-view-overlay-start sage-view-overlay-final nil nil nil))
      (overlay-put sage-view-current-overlay 'help-echo "Overlay made by sage-view")

      (let* ((base (expand-file-name (make-temp-name "output_")
				     (sage-view-temp-dir)))
	     (file (concat base ".tex")))
	(with-temp-file file
	  (insert sage-view-head)
	  (insert sage-view-insert)
	  (insert sage-view-tail))
	(sage-view-latex->pdf file)))))

(defun sage-view-output-filter-process-inline-plots (string)
  "Generate and place one overlay image for one inline plot,
found by looking for a particular png file in
`sage-view-temp-dir'.

This function expects the buffer to be narrowed to just the
current output; see `sage-view-output-filter' for how to do
that."

  ;; we need to foil emacs' image cache, which doesn't reload files with the same name
  (let* ((pngname (format "%s/sage-view.png" (sage-view-temp-dir)))
	 (base (expand-file-name (make-temp-name "plot_") (sage-view-temp-dir)))
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

      (goto-char (point-max))
      (let ((im (create-image pngname2 'png nil)))
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
      (narrow-to-region comint-last-input-end (process-mark (get-buffer-process (current-buffer))))

      (when sage-view-inline-plots-enabled ;; do we even want plot images inline?
	(sage-view-output-filter-process-inline-plots string))

      (when sage-view-inline-output-enabled ;; do we even want output latexed inline?
	(sage-view-output-filter-process-inline-output string)))))

(defun sage-view-gs-open ()
  "Start a Ghostscript conversion pass.")

(defun sage-view-temp-dir ()
  (unless (and sage-view-temp-dir-name
	       (file-exists-p sage-view-temp-dir-name)
	       (file-readable-p sage-view-temp-dir-name))
    (setq sage-view-temp-dir-name (make-temp-file (expand-file-name "tmp" "~/.sage/temp/") t))
    (unless (file-exists-p sage-view-temp-dir-name)
      (make-directory sage-view-temp-dir-name)))
  sage-view-temp-dir-name)

(defun sage-view-enable-inline-output ()
  (interactive)
  (setq sage-view-inline-output-enabled t)
  (python-send-receive-multiline "pretty_print_default(True);"))

(defun sage-view-disable-inline-output ()
  (interactive)
  (setq sage-view-inline-output-enabled nil)
  (python-send-receive-multiline "pretty_print_default(False);"))

(defun sage-view-enable-inline-plots ()
  (interactive)
  (setq sage-view-inline-plots-enabled t)
  (python-send-receive-multiline "sage.plot.plot.DOCTEST_MODE = True;")
  (python-send-receive-multiline (format "sage.plot.plot.DOCTEST_MODE_FILE = '%s';" (format "%s/sage-view.png" (sage-view-temp-dir)))))

(defun sage-view-disable-inline-plots ()
  (interactive)
  (setq sage-view-inline-plots-enabled nil)
  (python-send-receive-multiline "sage.plot.plot.DOCTEST_MODE = False;")
  (python-send-receive-multiline (format "sage.plot.plot.DOCTEST_MODE_FILE = None;")))

(define-minor-mode sage-view
  "With this mode, output in SAGE interactive buffers is
  preprocessed and texify." nil
  :group 'sage
  :lighter " Sage-View"
  (if sage-view
      (progn
	;; (sage-view-pretty-print-enable)
	(sage-view-enable-inline-output)
	(sage-view-enable-inline-plots)
	(make-local-variable 'comint-preoutput-filter-functions)
	(make-local-variable 'comint-output-filter-function)
	(add-hook 'comint-output-filter-functions 'sage-view-output-filter))
    (progn
      (sage-view-disable-inline-output)
      (sage-view-disable-inline-plots)
      (remove-hook 'comint-output-filter-functions 'sage-view-output-filter)
      (remove-hook 'comint-preoutput-filter-functions 'sage-view-preoutput-filter)
      (dired-delete-file (sage-view-temp-dir) 'always)
      (sage-view-pretty-print-disable)
      )))

(defun sage-view-always ()
  (sage-view 1))

(provide 'sage-view)
;;; sage-view.el ends here
