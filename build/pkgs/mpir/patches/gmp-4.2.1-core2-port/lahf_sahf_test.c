/*
 * Simple C code to see if the "lahf" and "sahf" instructions are
 * available.  The earliest amd64 processors don't have these
 * instructions, so we need to test for them to keep our code
 * backwards compatable for these older CPUs.
 */
#include<stdio.h>

/*
 * FUNCTION: int are_available_lahf_sahf(void)
 * Returns 1 if lahf and sahf instructions are available.  0 otherwise.
 */
int are_available_lahf_sahf(void)
{
  register unsigned long long rcx asm ("rcx");
  rcx = 0;
  asm ("movl\t$0x80000001,%eax");
  asm ("cpuid");
  return(rcx & 0x0000000000000001LL);
}

main()
{
  if (are_available_lahf_sahf())
    {
      printf("Yes");
    }
  else
    {
      printf("No");
    }
}
