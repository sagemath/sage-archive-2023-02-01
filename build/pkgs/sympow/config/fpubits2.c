/*
 * Copyright Jeroen Demeyer 2011
 *
 * Compute a*b + c
 * If the processor has a fused multiply-and-add instruction
 * (as on ia64), such an instruction will normally be used here.  Since
 * a fused multiply-and-add only rounds after the addition, it will
 * cause a larger apparent precision.
 *
 * We put this function in a separate file to make sure the compiler
 * does not optimize away the call to this function.
 */
double mul_and_add(double a, double b, double c)
{
        return a*b + c;
}
