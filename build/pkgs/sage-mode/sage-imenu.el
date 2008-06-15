;;; sage-test.el --- imenu for Sage code

;; Copyright (C) 2008  Nicholas Alexander

;; Author: Nicholas Alexander <ncalexan@pv109055.reshsg.uci.edu>
;; Keywords: sage imenu

(defun sage-imenu-create-index ())
;;   (let ((alist (python-imenu-create-index)))
;;     (

(defun sage-imenu-hook ()
  (setq imenu-create-index-function #'sage-imenu-create-index))

(add-hook 'sage-mode-hook 'sage-imenu-hook)

(provide 'sage-imenu)