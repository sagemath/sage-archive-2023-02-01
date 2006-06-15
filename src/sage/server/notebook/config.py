"""
Customization of the Notebook
"""

from js import JSKeyCode

JSKeyCode('request_completions', key = 39,  ctrl=True)
JSKeyCode('request_history',     key = 112, ctrl=True)
JSKeyCode('request_history_full',key = 114, ctrl=True)
JSKeyCode('close_helper',        key = 27)
JSKeyCode('interrupt',           key = 27)
JSKeyCode('send_input',          key = 13,  shift=True)
JSKeyCode('send_input_timed',    key = 13,  ctrl =True)
#JSKeyCode('prev_cell',           key = 38,  ctrl =True)
#JSKeyCode('next_cell',           key = 40,  ctrl =True)
JSKeyCode('delete_cell',         key = 8)
