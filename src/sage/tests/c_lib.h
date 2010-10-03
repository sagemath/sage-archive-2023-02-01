/*
 * C functions for use in sage/tests
 */

#ifndef SAGE_TESTS_C_LIB_H
#define SAGE_TESTS_C_LIB_H

#include <sys/types.h>

/* Wait ``ms`` milliseconds */
void ms_sleep(long ms);

/*
 * Wait ``ms`` milliseconds, then signal ``killpid`` with signal
 * ``signum``.  Wait ``interval`` milliseconds, then signal again.
 * Repeat this until ``n`` signals have been sent.  Usually, ``n``
 * will be equal to 1.  In that case, ``interval`` is irrelevant.
 */
void signal_pid_after_delay(int signum, pid_t killpid, long ms, long interval, int n);

/* Signal the Sage process */
#define signal_after_delay(signum, ms) signal_pid_after_delay(signum, getpid(), ms, 0, 1)

/* The same as above, but sending ``n`` signals */
#define signals_after_delay(signum, ms, interval, n) signal_pid_after_delay(signum, getpid(), ms, interval, n)


#endif
