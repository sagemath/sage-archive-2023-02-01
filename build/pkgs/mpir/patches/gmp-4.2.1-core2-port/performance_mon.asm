dnl
dnl File: performance_mon.s
dnl
dnl Description: A set of routines for monitoring timing, cache use,
dnl branch (mis)prediction, etc. of my code.
dnl

dnl
dnl This is just a check to see if we are in my code testing sandbox
dnl or if we are actually in the GMP source tree
dnl
ifdef(`__JWM_Test_Code__',`
include(`./config.m4')
define(`MPN_PREFIX',`jwm_mpn_')',`
include(`../config.m4')')


	.text
.globl GSYM_PREFIX()pm_read_time_stamp_counter
GSYM_PREFIX()pm_read_time_stamp_counter:
	rdtsc
	shl	$32,%rdx
	or	%rdx,%rax
	ret
