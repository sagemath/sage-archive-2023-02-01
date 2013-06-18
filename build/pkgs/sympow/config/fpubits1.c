/*
 * Copyright Jeroen Demeyer 2011
 *
 * Find out the precision of "double" floating-point numbers.
 * Return 0 if the precision is exactly 53,
 * return 1 if the precision is different from 53.
 */

#include <stdio.h>

double mul_and_add(double, double, double);  /* config/fpubits2.c */
void fpu_53bits();  /* src/fpu.c */

int main(int argc, char** argv)
{
	/* If x86 is defined, set the FPU to 53 bits */
#ifdef x86
	fpu_53bits();
#endif

        /* Let a = 1 + x,
         *     b = 1 + y,
         * and compute
         * a*b - 1 = (1 + x + y + x*y) - 1
         *         = x + y + x*y
         * The last equality will hold computationally if (1 + x + y + x*y)
         * can be represented exactly as floating-point number.
         * Assuming x and y are negative powers of 2, this will work as long
         * as x*y >= (1/2)^(p-1) where p is the precision of the FPU in bits.
         */

        int n = 0;
        double x = 1;
        double y = 1;
        for (;;)
        {
		/* Make sure x*y = (1/2)^n
		 * and that x and y are roughly equal. */
                n++;
                if (x >= y)
                        x /= 2;
                else
                        y /= 2;

                double r = mul_and_add(1+x, 1+y, -1);
                if (r != x + y + x*y) break;

                if (n >= 600)
                {
                        printf("The double precision of your FPU seems to be more than %i bits, bailing out.\n", n);
                        return 1;
                }
        }

        printf("The double precision of your FPU is %i bits.\n", n);
	if (n == 53) return 0;
	return 1;
}
