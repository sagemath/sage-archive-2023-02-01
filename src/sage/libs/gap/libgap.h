#include <gap/gapstate.h>
#include <gap/system.h>


#define _GAP_Error_Setjmp() !sySetjmp(STATE(ReadJmpError))
