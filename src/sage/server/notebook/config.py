"""
Customization of the Notebook
"""

#from js import JSKeyHandler
import js


js.keyhandler.add('request_introspections', key = 32,  ctrl=True)  # control right arrow
js.keyhandler.add('request_introspections', key = 9)  # tab
#js.keyhandler.add('request_introspections', key = 112, ctrl=True)  # F1
js.keyhandler.add('request_history',     key = 112, ctrl=True)
js.keyhandler.add('request_log',         key = 113, ctrl=True)
#js.keyhandler.add('request_history_full',key = 114, ctrl=True)
js.keyhandler.add('close_helper',        key = 27)
js.keyhandler.add('interrupt',           key = 27)
js.keyhandler.add('send_input',          key = 13,  shift=True)
js.keyhandler.add('send_input_newcell',  key = 13,  ctrl =True)
js.keyhandler.add('prev_cell',           key = 38,  ctrl =True)
js.keyhandler.add('next_cell',           key = 40,  ctrl =True)
js.keyhandler.add('page_up',           key = 33)
js.keyhandler.add('page_down',           key = 34)
js.keyhandler.add('delete_cell',         key = 8)
js.keyhandler.add('generic_submit',      key = 13)

js.keyhandler.add('menu_left',           key = 37)
js.keyhandler.add('menu_up',             key = 38)
js.keyhandler.add('menu_right',          key = 39)
js.keyhandler.add('menu_down',           key = 40)
js.keyhandler.add('menu_pick',           key = 13)


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
