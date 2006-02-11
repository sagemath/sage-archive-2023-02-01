#include "interrupt.h"

void sig_handle(int n) {
  signal(SIGINT, last_handler);  /* restore python signal handler (or will crash when
                                    user presses ctrl-c in the interpreter. */
  siglongjmp(env, n);
}
