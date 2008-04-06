"""nodoctest
Customization of the Notebook
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
js.keyhandler.add('enter',                  key = "KEY_ENTER")
js.keyhandler.add('enter',                  key = "KEY_RETURN")
js.keyhandler.add('split_cell',             key = "KEY_ENTER", ctrl=True)
js.keyhandler.add('split_cell',             key = "KEY_RETURN",ctrl=True)
js.keyhandler.add('split_cell',             key = "KEY_CTRLENTER", ctrl=True)  # needed on OS X Firefox
js.keyhandler.add('join_cell',              key = "KEY_BKSPC", ctrl=True)

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
