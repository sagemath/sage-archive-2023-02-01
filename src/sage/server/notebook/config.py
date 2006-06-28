"""
Customization of the Notebook
"""

from js import JSKeyCode

#JSKeyCode('request_introspections', key = 39,  ctrl=True)  # control right arrow
JSKeyCode('request_introspections', key = 9)  # tab
#JSKeyCode('request_introspections', key = 112)  # F1
JSKeyCode('request_history',     key = 112, ctrl=True)
JSKeyCode('request_log',         key = 113, ctrl=True)
#JSKeyCode('request_history_full',key = 114, ctrl=True)
JSKeyCode('close_helper',        key = 27)
JSKeyCode('interrupt',           key = 27)
JSKeyCode('send_input',          key = 13,  shift=True)
JSKeyCode('send_input_newcell',  key = 13,  ctrl =True)
JSKeyCode('prev_cell',           key = 38,  ctrl =True)
JSKeyCode('next_cell',           key = 40,  ctrl =True)
JSKeyCode('page_up',           key = 33)
JSKeyCode('page_down',           key = 34)
JSKeyCode('delete_cell',         key = 8)
JSKeyCode('generic_submit',      key = 13)
JSKeyCode('menu_prev',           key = 38)
JSKeyCode('menu_next',           key = 40)
JSKeyCode('menu_pick',           key = 13)


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
