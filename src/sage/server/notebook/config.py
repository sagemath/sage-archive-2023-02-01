"""
Notebook Keybindings

This module is responsible for setting the keyboard bindings for the notebook.

These are the standard key and mouse bindings available in the
notebook:

- *Evaluate Input:* Press **shift-enter**. You can start several calculations at once. If you press **alt-enter** instead, then a new cell is created after the current one. If you press **ctrl-enter** then the cell is split and both pieces are evaluated separately.

- *Tab Completion:* Press **tab** while the cursor is on an identifier. On some web browsers (e.g., Opera) you must use control-space instead of tab.

- *Insert New Cell:* Put the mouse between an output and input until the horizontal line appears and click. If you press Alt-Enter in a cell, the cell is evaluated and a new cell is inserted after it.

- *Delete Cell:* Delete all cell contents, then press **backspace**.

- *Split and Join Cells:* Press **ctrl-;** in a cell to split it into two cells, and **ctrl-backspace** to join them. Press **ctrl-enter** to split a cell and evaluate both pieces.

- *Insert New HTML Cell:* Shift click between cells to create a new HTML cell. Double click on existing HTML to edit it. Use $...$ and $$...$$ to include typeset math in the HTML block.

- *Hide/Show Output:* Click on the left side of output to toggle between hidden, shown with word wrap, and shown without word wrap.

- *Indenting Blocks:* Highlight text and press **>** to indent it all and **<** to unindent it all (works in Safari and Firefox). In Firefox you can also press tab and shift-tab.

- *Comment/Uncomment Blocks:* Highlight text and press **ctrl-.** to comment it and **ctrl-,** to uncomment it. Alternatively, use **ctrl-3** and **ctrl-4**.

"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

#from js import JSKeyHandler
import js


js.keyhandler.add('request_introspections', key = "KEY_SPC",   ctrl=True)  # control space
js.keyhandler.add('request_introspections', key = "KEY_TAB",   shift=False)  # tab
js.keyhandler.add('indent',                 key = "KEY_TAB",   shift=False)  # tab
js.keyhandler.add('indent',                 key = "KEY_GT",    shift=False)  # tab
js.keyhandler.add('unindent',               key = "KEY_TAB",   shift=True)  # tab
js.keyhandler.add('unindent',               key = "KEY_LT",    shift=False)
js.keyhandler.add('request_history',        key = "KEY_Q",     ctrl=True)
js.keyhandler.add('request_history',        key = "KEY_QQ",    ctrl=True)
js.keyhandler.add('request_log',            key = "KEY_P",     ctrl=True)
js.keyhandler.add('request_log',            key = "KEY_PP",    ctrl=True)
js.keyhandler.add('close_helper',           key = "KEY_ESC")
js.keyhandler.add('interrupt',              key = "KEY_ESC")
js.keyhandler.add('send_input',             key = "KEY_RETURN", shift=True)
js.keyhandler.add('send_input',             key = "KEY_ENTER", shift=True)
js.keyhandler.add('send_input_newcell',     key = "KEY_RETURN", alt =True)
js.keyhandler.add('send_input_newcell',     key = "KEY_ENTER", alt =True)
js.keyhandler.add('prev_cell',              key = "KEY_UP",    ctrl =True)
js.keyhandler.add('next_cell',              key = "KEY_DOWN",  ctrl =True)
js.keyhandler.add('page_up',                key = "KEY_PGUP")
js.keyhandler.add('page_down',              key = "KEY_PGDN")
js.keyhandler.add('delete_cell',            key = "KEY_BKSPC")
js.keyhandler.add('generic_submit',         key = "KEY_ENTER")
js.keyhandler.add('up_arrow',               key = "KEY_UP")
js.keyhandler.add('down_arrow',             key = "KEY_DOWN")
js.keyhandler.add('comment',                key = "KEY_DOT",   ctrl=True)
js.keyhandler.add('uncomment',              key = "KEY_COMMA", ctrl=True)
js.keyhandler.add('comment',                key = "KEY_3",     ctrl=True)
js.keyhandler.add('uncomment',              key = "KEY_4",     ctrl=True)

js.keyhandler.add('control',                key = "KEY_CTRL")
js.keyhandler.add('backspace',              key = "KEY_BKSPC")
js.keyhandler.add('enter',                  key = "KEY_ENTER")
js.keyhandler.add('enter',                  key = "KEY_RETURN")
js.keyhandler.add('enter_shift',            key = "KEY_ENTER", shift=True)
js.keyhandler.add('enter_shift',            key = "KEY_RETURN", shift=True)
js.keyhandler.add('spliteval_cell',         key = "KEY_ENTER", ctrl=True)
js.keyhandler.add('spliteval_cell',         key = "KEY_RETURN",ctrl=True)
js.keyhandler.add('spliteval_cell',         key = "KEY_CTRLENTER", ctrl=True)  # needed on OS X Firefox
js.keyhandler.add('join_cell',              key = "KEY_BKSPC", ctrl=True)
js.keyhandler.add('split_cell',             key = "KEY_SEMI", ctrl=True)
js.keyhandler.add('split_cell_noctrl',      key = "KEY_SEMI")

js.keyhandler.add('menu_left',           key = "KEY_LEFT")
js.keyhandler.add('menu_up',             key = "KEY_UP")
js.keyhandler.add('menu_right',          key = "KEY_RIGHT")
js.keyhandler.add('menu_down',           key = "KEY_DOWN")
js.keyhandler.add('menu_pick',           key = "KEY_ENTER")
js.keyhandler.add('menu_pick',           key = "KEY_RETURN")


"""
8  -- backspace
9  -- tab
13 -- return
27 -- esc
32 -- space
37 -- left
38 -- up
39 -- right
40 -- down
"""
