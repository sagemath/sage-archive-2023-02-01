#include <pari/pari.h>
#include "interrupt.h"


/*****************************************
   Interrupts and PARI exception handling
 *****************************************/

#define pari_catch_sig_on() _sig_on_(NULL)
#define pari_catch_sig_str(s) _sig_on_(s)
#define pari_catch_sig_off() (_sig_off_(__FILE__, __LINE__), _pari_check_warning())
