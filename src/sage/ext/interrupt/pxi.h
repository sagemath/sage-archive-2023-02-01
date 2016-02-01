/*
 * The order of these includes is very important, as each uses
 * stuff defined by the previous. That is also the reason why
 * struct_signals.h and macros.h must be separate files.
 */
#include "interrupt/struct_signals.h"
#include "interrupt/interrupt_api.h"
#include "interrupt/macros.h"

/* Undefine this macro from interrupt_api.h to avoid compiler warnings:
 * Cython redefines it when cimporting interrupt.pxd */
#undef _signals
