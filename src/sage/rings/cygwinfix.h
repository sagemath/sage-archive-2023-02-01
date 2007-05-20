// cygwinfix.h
// defines NAN and INFINITY on cygwin

#if defined(__CYGWIN__)
#define NAN (0.0/0.0)
#define INFINITY HUGE_VALF
#endif
