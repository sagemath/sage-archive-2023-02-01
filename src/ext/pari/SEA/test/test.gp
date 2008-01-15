read("sea");
allocatemem(200*10^6);
{
  min_size = 10;
  max_size = 300;
  nbr_tests = 100;
  timings = vector(max_size);
  for(i = 1, nbr_tests,
    forstep (nb = min_size, max_size, 10,
      p = nextprime(random(2^nb)); print("p = ", p);
      print("i = ",i);
      A = random(p); print("A = ", A);
      B = random(p); print("B = ", B);
      if (is_singular(A, B, p), next);
      E = ellinit([0, 0, 0, A, B]*Mod(1,p));
      t_i = gettime(); ellsea(E,p);
      t_f = gettime();
      timings[nb] = round(((i-1)*timings[nb]+t_f-t_i+0.)/i)
    )
  )
}
