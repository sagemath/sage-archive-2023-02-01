


/////////////////////////////////////////////////
// Computes the absolute value of an mpq_class //
/////////////////////////////////////////////////

mpq_class Abs(mpq_class m) {
  if (m > 0)
    return m;
  else
    return -m;
}





////////////////////////////////////////
// Raise a rational number to a power //
////////////////////////////////////////

mpq_class RationalPower(mpq_class base, unsigned int pow) {

  mpz_t top, bot;
  mpz_init (top);
  mpz_init (bot);

  mpq_get_num(top, base.get_mpq_t());
  mpq_get_den(bot, base.get_mpq_t());

  mpz_pow_ui(top, top, pow);
  mpz_pow_ui(bot, bot, pow);


  mpz_class top1, bot1;
  top1 = mpz_class(top);
  bot1 = mpz_class(bot);

  // Unallocate the mpz_t variables
  mpz_clear(top);
  mpz_clear(bot);


  return mpq_class(top1, bot1);
}


////////////////////////////////////////////////////////////////////////////
// Take the square root of a rational number already known to be a square //
////////////////////////////////////////////////////////////////////////////

mpq_class RationalSqrt(mpq_class m) {

  mpz_t top, bot;
  mpz_init (top);
  mpz_init (bot);

  mpq_get_num(top, m.get_mpq_t());
  mpq_get_den(bot, m.get_mpq_t());

  mpz_sqrt(top, top);
  mpz_sqrt(bot, bot);


  mpz_class top1, bot1;
  top1 = mpz_class(top);
  bot1 = mpz_class(bot);

  // Unallocate the mpz_t variables
  mpz_clear(top);
  mpz_clear(bot);


  return mpq_class(top1, bot1);
}

