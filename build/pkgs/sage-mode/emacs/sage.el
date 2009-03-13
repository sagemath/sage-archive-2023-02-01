;;;_* sage.el --- Major modes for editing sage code, interacting with inferior
;;;sage processes, building sage, and doctesting sage

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

;; sage.el provides customization options and autoloads.

;;; Code:

(defcustom inferior-sage-prompt (rx line-start (1+ (and (or "sage:" "....." ">>>" "..." "(Pdb)" "ipdb>" "(gdb)") " ")))
  "Regular expression matching the SAGE prompt."
  :group 'sage
  :type 'regexp)

(defcustom inferior-sage-timeout 500000
  "How long to wait for a SAGE prompt."
  :group 'sage
  :type 'integer)

(defcustom sage-command (expand-file-name "~/bin/sage")
  "Actual command used to run SAGE.
Additional arguments are added when the command is used by `run-sage' et al."
  :group 'sage
  :type 'string)

(defcustom sage-startup-command "import sage_emacs as emacs"
  "Run this command each time SAGE slave is executed by `run-sage'."
  :group 'sage
  :type 'string)

(defcustom sage-startup-hook nil
  "*Normal hook (list of functions) run after `sage' is run and the first prompt is seen.
See `run-hooks'."
  :group 'sage
  :type 'hook)

(defgroup sage-test nil "Run Sage doctests"
  :group 'sage
  :prefix "sage-test-"
  :link '(url-link :tag "Homepage" "http://wiki.sagemath.org/sage-mode"))

(defcustom sage-test-setup-hook nil
  "List of hook functions run by `sage-test-process-setup' (see `run-hooks')."
  :type 'hook
  :group 'sage-test)

(defcustom sage-test-prompt (rx line-start (? (or "-" "+")) (0+ (or space punct)) (1+ (or "sage: " ">>> " "....." "...")))
  "Regular expression matching the SAGE prompt of a single doctest line."
  :group 'sage-test
  :type 'regexp)

(defgroup sage-build nil "Build the Sage library"
  :group 'sage
  :prefix "sage-build-"
  :link '(url-link :tag "Homepage" "http://wiki.sagemath.org/sage-mode"))

(defcustom sage-build-setup-hook nil
  "List of hook functions run by `sage-build-process-setup' (see `run-hooks')."
  :type 'hook
  :group 'sage-build)

(defcustom sage-rerun-command (format "%s" sage-command)
  "Actual command used to rerun SAGE.
Additional arguments are added when the command is used by `rerun-sage' et al."
  :group 'sage-build
  :type 'string)

(defgroup sage-view nil "Typeset Sage output on the fly"
  :group 'sage
  :prefix "sage-view-"
  :link '(url-link :tag "Homepage" "http://wiki.sagemath.org/sage-mode"))

(defcustom sage-view-gs-command (if (eq system-type 'windows-nt)
				    "GSWIN32C.EXE"
				  "gs")
  "*Ghostscript command to convert from PDF to PNG.

See also `sage-view-gs-options', `sage-view-anti-aliasing-level'
and `sage-view-scale'."
  :type 'string
  :group 'sage-view)

(defcustom sage-view-gs-options
  '("-sDEVICE=png16m" "-dBATCH" "-dSAFER" "-q" "-dNOPAUSE")
  "*Options for Ghostscript when converting from PDF to PNG."
  :type 'list
  :group 'sage-view)

(defcustom sage-view-anti-aliasing-level 2
  "*Level of anti-aliasing used when converting from PDF to PNG. "
  :type 'integer
  :group 'sage-view)

(defcustom sage-view-scale 1.2
  "*Scale used when converting from PDF/PS to PNG."
  :type 'number
  :group 'sage-view)

(defcustom sage-view-margin '(1 . 1)
  "*Margin (in pixels or (pixels-x . pixels-y)) added around displayed images."
  :type '(choice integer (cons integer integer))
  :group 'sage-view)

(defcustom sage-view-scale-factor 0.2
  "*Factor used when zooming."
  :type 'number
  :group 'sage-view)

;; Generate autoload for sage elisp files
(defun sage-update-autoloads nil
  "Generate autoloads for sage elisp files.

WARNING: do not use this unless you are distributing a new
version of `sage-mode'!"
  (interactive)
  (let ((generated-autoload-file "~/emacs/sage/emacs/sage-load.el"))
    (update-directory-autoloads "~/emacs/sage/emacs")))

(load "sage-load")
(provide 'sage)
