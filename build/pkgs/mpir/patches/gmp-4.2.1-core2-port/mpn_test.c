/* **************************************************************************
   File: mpn_test.c

************************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define NUM_TRIALS 500
#define MAX_NUM_LIMBS 5000
#define SANE_LIMB_TEST 1248
#define RANDOM_SEED 1234
#define RANDOM_SEED2 93857263
#define RANDOM_SEED3 592023


/* modlimb_invert_table[i] is the multiplicative inverse of 2*i+1 mod 256,
   ie. (modlimb_invert_table[i] * (2*i+1)) % 256 == 1 */

const unsigned char  modlimb_invert_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
};




uint64_t jwm_mpn_add_n(uint64_t*,uint64_t*,uint64_t*,uint64_t);
uint64_t orig_mpn_add_n(uint64_t*,uint64_t*,uint64_t*,uint64_t);

uint64_t jwm_mpn_sub_n(uint64_t*,uint64_t*,uint64_t*,uint64_t);
uint64_t orig_mpn_sub_n(uint64_t*,uint64_t*,uint64_t*,uint64_t);

uint64_t jwm_mpn_addmul_1(uint64_t*,uint64_t*,uint64_t,uint64_t);
uint64_t orig_mpn_addmul_1(uint64_t*,uint64_t*,uint64_t,uint64_t);

uint64_t pm_read_time_stamp_counter();

int compare_uint64_t(const void *p1,const void *p2)
{
  uint64_t value1,value2;
  uint64_t difference;
  int ecode;

  value1 = *((uint64_t *)p1);
  value2 = *((uint64_t *)p2);
  ecode = (int)(value1 - value2);

  return(ecode);
}

int compare_big_ints(uint64_t *p1,uint64_t *p2,uint64_t num_limbs)
{
  uint64_t i;
  int ecode = 1;

  for (i=0;i<num_limbs;i++)
    {
      if (compare_uint64_t(p1+i,p2+i) != 0)
	{
	  ecode = 0;
	  break;
	}
    }
  return(ecode);
}

void random_init_big_int(uint32_t seed, uint64_t *p, uint64_t num_limbs)
{
  uint32_t randha,randhb,randla,randlb,randh,randl;
  uint64_t value,i;

  srandom(seed);

  for (i=0;i<num_limbs;i++)
    {
      randla = random();
      randlb = random();
      randha = random();
      randhb = random();
      randh = (randha << 4) ^ randhb;
      randl = (randla << 4) ^ randlb;
      value = ( ((uint64_t)randh) << 32) | ((uint64_t)randl);
      p[i] = value;
    }
  return;
}





int sub_n_test(uint64_t* results,uint64_t num_limbs)
{
  uint64_t big_rand_int1[MAX_NUM_LIMBS];
  uint64_t big_rand_int2[MAX_NUM_LIMBS];
  uint64_t big_result_my[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_old[2*MAX_NUM_LIMBS+1];
  uint64_t orig_mpn_trial_times[NUM_TRIALS];
  uint64_t jwm_mpn_trial_times[NUM_TRIALS];
  uint64_t carry_out_my, carry_out_old, e_code;
  uint64_t time_stamp1, time_stamp2;
  int i;

  //
  // Initialize our big_rand_ints
  //
  random_init_big_int(RANDOM_SEED,big_rand_int1,num_limbs);
  random_init_big_int(RANDOM_SEED2,big_rand_int2,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);


  //
  // First a correctness test:
  //

  carry_out_old = orig_mpn_sub_n(big_result_old,
				 big_rand_int1,
				 big_rand_int2,
				 num_limbs);


  carry_out_my = jwm_mpn_sub_n(big_result_my,
			       big_rand_int1,
			       big_rand_int2,
			       num_limbs);

  if ( !(compare_big_ints(big_result_old,
			  big_result_my,
			  num_limbs) && (carry_out_old == carry_out_my) ) )
    {
      for(i=0;i<num_limbs;i++)
	{
	  printf("big_rand_int1: %.16llx   big_rand_int2: %.16llx\n",
		 big_rand_int1[i],
		 big_rand_int2[i]);
	}

      for(i=0;i<num_limbs;i++)
	{
	  printf("jwm_mpn_sub_n: %.16llx      orig_mpn_sub_n: %.16llx\n",
		 big_result_my[i],
		 big_result_old[i]);
	}
      printf("my carry: %.16llx              old carry: %.16llx\n",
	     carry_out_my,
	     carry_out_old);
      printf("ERROR -- adds do not match!\n\n");
      exit(-1);
    }

  //
  // Now for timing tests
  //

  for(i=0;i<NUM_TRIALS;i++)
    {
      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_my = jwm_mpn_sub_n(big_result_my,
				   big_rand_int1,
				   big_rand_int2,
				   num_limbs);

      time_stamp2 = pm_read_time_stamp_counter();
      jwm_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(jwm_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  for(i=0;i<NUM_TRIALS;i++)
    {
      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_old = orig_mpn_sub_n(big_result_old,
				     big_rand_int1,
				     big_rand_int2,
				     num_limbs);

      time_stamp2 = pm_read_time_stamp_counter();
      orig_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(orig_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  results[0] = jwm_mpn_trial_times[0];
  results[1] = orig_mpn_trial_times[0];
  return(1);
}



int add_n_test(uint64_t* results,uint64_t num_limbs)
{
  uint64_t big_rand_int1[MAX_NUM_LIMBS];
  uint64_t big_rand_int2[MAX_NUM_LIMBS];
  uint64_t big_result_my[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_old[2*MAX_NUM_LIMBS+1];
  uint64_t orig_mpn_trial_times[NUM_TRIALS];
  uint64_t jwm_mpn_trial_times[NUM_TRIALS];
  uint64_t carry_out_my, carry_out_old, e_code;
  uint64_t time_stamp1, time_stamp2;
  int i;

  //
  // Initialize our big_rand_ints
  //
  random_init_big_int(RANDOM_SEED,big_rand_int1,num_limbs);
  random_init_big_int(RANDOM_SEED2,big_rand_int2,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);


  //
  // First a correctness test:
  //

  carry_out_old = orig_mpn_add_n(big_result_old,
				 big_rand_int1,
				 big_rand_int2,
				 num_limbs);


  carry_out_my = jwm_mpn_add_n(big_result_my,
			       big_rand_int1,
			       big_rand_int2,
			       num_limbs);

  if ( !(compare_big_ints(big_result_old,
			  big_result_my,
			  num_limbs) && (carry_out_old == carry_out_my) ) )
    {
      for(i=0;i<num_limbs;i++)
	{
	  printf("big_rand_int1: %.16llx   big_rand_int2: %.16llx\n",
		 big_rand_int1[i],
		 big_rand_int2[i]);
	}

      for(i=0;i<num_limbs;i++)
	{
	  printf("jwm_mpn_add_n: %.16llx      orig_mpn_add_n: %.16llx\n",
		 big_result_my[i],
		 big_result_old[i]);
	}
      printf("my carry: %.16llx              old carry: %.16llx\n",
	     carry_out_my,
	     carry_out_old);
      printf("ERROR -- adds do not match!\n\n");
      exit(-1);
    }

  //
  // Now for timing tests
  //

  for(i=0;i<NUM_TRIALS;i++)
    {
      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_my = jwm_mpn_add_n(big_result_my,
				   big_rand_int1,
				   big_rand_int2,
				   num_limbs);

      time_stamp2 = pm_read_time_stamp_counter();
      jwm_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(jwm_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  for(i=0;i<NUM_TRIALS;i++)
    {
      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_old = orig_mpn_add_n(big_result_old,
				     big_rand_int1,
				     big_rand_int2,
				     num_limbs);

      time_stamp2 = pm_read_time_stamp_counter();
      orig_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(orig_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  results[0] = jwm_mpn_trial_times[0];
  results[1] = orig_mpn_trial_times[0];
  return(1);
}



int addmul_1_test(uint64_t* results,uint64_t num_limbs)
{
  uint64_t big_rand_int1[MAX_NUM_LIMBS];
  uint64_t big_rand_int2[MAX_NUM_LIMBS];
  uint64_t big_result_my[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_my_before[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_old[2*MAX_NUM_LIMBS+1];
  uint64_t orig_mpn_trial_times[NUM_TRIALS];
  uint64_t jwm_mpn_trial_times[NUM_TRIALS];
  uint64_t carry_out_my, carry_out_old, e_code;
  uint64_t time_stamp1, time_stamp2;
  int i;

  //
  // Initialize our big_rand_ints
  //
  random_init_big_int(RANDOM_SEED,big_rand_int1,num_limbs);
  random_init_big_int(RANDOM_SEED2,big_rand_int2,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);


  //
  // First a correctness test:
  //


  carry_out_my = jwm_mpn_addmul_1(big_result_my,
				  big_rand_int1,
				  num_limbs,
				  big_rand_int2[num_limbs-1]);

  carry_out_old = orig_mpn_addmul_1(big_result_old,
				    big_rand_int1,
				    num_limbs,
				    big_rand_int2[num_limbs-1]);


  if ( !(compare_big_ints(big_result_old,
			  big_result_my,
			  num_limbs) && (carry_out_my == carry_out_old) ) )
    {
      random_init_big_int(RANDOM_SEED3,big_result_my_before,num_limbs);

      printf("ERROR ---- ERROR ---- ERROR\nInputs:\n");
      printf("s2limb:\t\t%.16llx\n\n",big_rand_int2[num_limbs-1]);
      for(i=0;i<num_limbs;i++)
	{
	  printf("big_rand_int1:\t\t%.16llx\tbig_result_before:\t\t%.16llx\n",
		 big_rand_int1[i],
		 big_result_my_before[i]);
	}
      printf("\nResults:\n");
      for(i=0;i<num_limbs;i++)
	{
	  printf("jwm_mpn_addmul_1:\t%.16llx\torig_mpn_addmul_1:\t%.16llx\n",
		 big_result_my[i],
		 big_result_old[i]);
	}
      printf("\nCarry Out:\n");
      printf("my carry:\t%.16llx\t\told carry:\t%.16llx\n",
	     carry_out_my,
	     carry_out_old);
      printf("ERROR -- addmul_1 does not match!\n\n");
      exit(-1);
    }

  //
  // Now for timing tests
  //

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_my = jwm_mpn_addmul_1(big_result_my,
				      big_rand_int1,
				      num_limbs,
				      big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      jwm_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(jwm_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_old = orig_mpn_addmul_1(big_result_old,
					big_rand_int1,
					num_limbs,
					big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      orig_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(orig_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  results[0] = jwm_mpn_trial_times[0];
  results[1] = orig_mpn_trial_times[0];
  return(1);
}



int submul_1_test(uint64_t* results,uint64_t num_limbs)
{
  uint64_t big_rand_int1[MAX_NUM_LIMBS];
  uint64_t big_rand_int2[MAX_NUM_LIMBS];
  uint64_t big_result_my[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_my_before[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_old[2*MAX_NUM_LIMBS+1];
  uint64_t orig_mpn_trial_times[NUM_TRIALS];
  uint64_t jwm_mpn_trial_times[NUM_TRIALS];
  uint64_t carry_out_my, carry_out_old, e_code;
  uint64_t time_stamp1, time_stamp2;
  int i;

  //
  // Initialize our big_rand_ints
  //
  random_init_big_int(RANDOM_SEED,big_rand_int1,num_limbs);
  random_init_big_int(RANDOM_SEED2,big_rand_int2,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);


  //
  // First a correctness test:
  //


  carry_out_my = jwm_mpn_submul_1(big_result_my,
				  big_rand_int1,
				  num_limbs,
				  big_rand_int2[num_limbs-1]);

  carry_out_old = orig_mpn_submul_1(big_result_old,
				    big_rand_int1,
				    num_limbs,
				    big_rand_int2[num_limbs-1]);


  if ( !(compare_big_ints(big_result_old,
			  big_result_my,
			  num_limbs) && (carry_out_my == carry_out_old) ) )
    {
      random_init_big_int(RANDOM_SEED3,big_result_my_before,num_limbs);

      printf("ERROR ---- ERROR ---- ERROR\nInputs:\n");
      printf("s2limb:\t\t%.16llx\n\n",big_rand_int2[num_limbs-1]);
      for(i=0;i<num_limbs;i++)
	{
	  printf("big_rand_int1:\t\t%.16llx\tbig_result_before:\t\t%.16llx\n",
		 big_rand_int1[i],
		 big_result_my_before[i]);
	}
      printf("\nResults:\n");
      for(i=0;i<num_limbs;i++)
	{
	  printf("jwm_mpn_submul_1:\t%.16llx\torig_mpn_submul_1:\t%.16llx\n",
		 big_result_my[i],
		 big_result_old[i]);
	}
      printf("\nCarry Out:\n");
      printf("my carry:\t%.16llx\t\told carry:\t%.16llx\n",
	     carry_out_my,
	     carry_out_old);
      printf("ERROR -- submul_1 does not match!\n\n");
      exit(-1);
    }

  //
  // Now for timing tests
  //

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_my = jwm_mpn_submul_1(big_result_my,
				      big_rand_int1,
				      num_limbs,
				      big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      jwm_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(jwm_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_old = orig_mpn_submul_1(big_result_old,
					big_rand_int1,
					num_limbs,
					big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      orig_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(orig_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  results[0] = jwm_mpn_trial_times[0];
  results[1] = orig_mpn_trial_times[0];
  return(1);
}



int mul_1_test(uint64_t* results,uint64_t num_limbs)
{
  uint64_t big_rand_int1[MAX_NUM_LIMBS];
  uint64_t big_rand_int2[MAX_NUM_LIMBS];
  uint64_t big_result_my[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_my_before[2*MAX_NUM_LIMBS+1];
  uint64_t big_result_old[2*MAX_NUM_LIMBS+1];
  uint64_t orig_mpn_trial_times[NUM_TRIALS];
  uint64_t jwm_mpn_trial_times[NUM_TRIALS];
  uint64_t carry_out_my, carry_out_old, e_code;
  uint64_t time_stamp1, time_stamp2;
  int i;

  //
  // Initialize our big_rand_ints
  //
  random_init_big_int(RANDOM_SEED,big_rand_int1,num_limbs);
  random_init_big_int(RANDOM_SEED2,big_rand_int2,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);
  random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);


  //
  // First a correctness test:
  //


  carry_out_my = jwm_mpn_mul_1(big_result_my,
				  big_rand_int1,
				  num_limbs,
				  big_rand_int2[num_limbs-1]);

  carry_out_old = orig_mpn_mul_1(big_result_old,
				    big_rand_int1,
				    num_limbs,
				    big_rand_int2[num_limbs-1]);


  if ( !(compare_big_ints(big_result_old,
			  big_result_my,
			  num_limbs) && (carry_out_my == carry_out_old) ) )
    {
      random_init_big_int(RANDOM_SEED3,big_result_my_before,num_limbs);

      printf("ERROR ---- ERROR ---- ERROR\nInputs:\n");
      printf("s2limb:\t\t%.16llx\n\n",big_rand_int2[num_limbs-1]);
      for(i=0;i<num_limbs;i++)
	{
	  printf("big_rand_int1:\t\t%.16llx\tbig_result_before:\t\t%.16llx\n",
		 big_rand_int1[i],
		 big_result_my_before[i]);
	}
      printf("\nResults:\n");
      for(i=0;i<num_limbs;i++)
	{
	  printf("jwm_mpn_mul_1:\t%.16llx\torig_mpn_mul_1:\t%.16llx\n",
		 big_result_my[i],
		 big_result_old[i]);
	}
      printf("\nCarry Out:\n");
      printf("my carry:\t%.16llx\t\told carry:\t%.16llx\n",
	     carry_out_my,
	     carry_out_old);
      printf("ERROR -- mul_1 does not match!\n\n");
      exit(-1);
    }

  //
  // Now for timing tests
  //

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_my,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_my = jwm_mpn_mul_1(big_result_my,
				      big_rand_int1,
				      num_limbs,
				      big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      jwm_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(jwm_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  for(i=0;i<NUM_TRIALS;i++)
    {
      random_init_big_int(RANDOM_SEED3,big_result_old,num_limbs);

      time_stamp1 = pm_read_time_stamp_counter();

      carry_out_old = orig_mpn_mul_1(big_result_old,
					big_rand_int1,
					num_limbs,
					big_rand_int2[num_limbs-1]);

      time_stamp2 = pm_read_time_stamp_counter();
      orig_mpn_trial_times[i] = time_stamp2-time_stamp1;
    }

  qsort(orig_mpn_trial_times,
	NUM_TRIALS,
	sizeof(uint64_t),
	compare_uint64_t);

  results[0] = jwm_mpn_trial_times[0];
  results[1] = orig_mpn_trial_times[0];
  return(1);
}



main()
{
  uint64_t results[2];
  uint64_t jwm_cycle_count, orig_cycle_count;
  uint64_t num_limbs;

  for(num_limbs=1;num_limbs<SANE_LIMB_TEST;num_limbs++)
    {
      add_n_test(results,num_limbs);
      printf("Limb count: \t\t%lli\n",num_limbs);
      printf("My add_n clock cycles: \t\t%lli \t clocks/limb: %f\n",
	     results[0],
	     ((double)results[0])/num_limbs);
      printf("Orig add_n clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[1],
	     ((double)results[1])/num_limbs);

      sub_n_test(results,num_limbs);
      printf("My sub_n clock cycles: \t\t%lli \t clocks/limb: %f\n",
	     results[0],
	     ((double)results[0])/num_limbs);
      printf("Orig sub_n clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[1],
	     ((double)results[1])/num_limbs);

      addmul_1_test(results,num_limbs);
      printf("My addmul_1 clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[0],
	     ((double)results[0])/num_limbs);
      printf("Orig addmul_1 clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[1],
	     ((double)results[1])/num_limbs);

      submul_1_test(results,num_limbs);
      printf("My submul_1 clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[0],
	     ((double)results[0])/num_limbs);
      printf("Orig submul_1 clock cycles: \t%lli \t clocks/limb: %f\n",
	     results[1],
	     ((double)results[1])/num_limbs);

      mul_1_test(results,num_limbs);
      printf("My mul_1 clock cycles: \t\t%lli \t clocks/limb: %f\n",
	     results[0],
	     ((double)results[0])/num_limbs);
      printf("Orig mul_1 clock cycles: \t%lli \t clocks/limb: %f\n\n",
	     results[1],
	     ((double)results[1])/num_limbs);

    }
  return(1);
}
