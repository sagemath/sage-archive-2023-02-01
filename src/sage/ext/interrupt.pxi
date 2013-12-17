#
# See c_lib/include/interrupt.h
#
cdef extern from 'interrupt.h':
    int sig_on() except 0
    int sig_str(char*) except 0
    int sig_check() except 0
    int sig_on_no_except()
    int sig_str_no_except(char*)
    int sig_check_no_except()
    void sig_off()
    void sig_retry()
    void sig_error()
    void sig_block()
    void sig_unblock()
    void set_sage_signal_handler_message(char* s)

    ctypedef struct sage_signals_t:
        int sig_on_count
        int interrupt_received
        int inside_signal_handler
        int block_sigint
        char* s
        int (*raise_exception)(int sig, const char* msg) except 0

    sage_signals_t _signals

    # These provide backwards compatibility with sage-4.6 and earlier
    int _sig_on
    void _sig_str(char*)
    int _sig_off

    void cython_check_exception() except *
