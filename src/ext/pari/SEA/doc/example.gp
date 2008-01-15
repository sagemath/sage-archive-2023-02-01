setrand(1);
read("sea");
ellrandom(bit) =
{ local(p,a4,a6);
  p = nextprime(random(1<<bit));
  a4 = random(p);
  a6 = random(p);
  ellinit([0,0,0,a4,a6]*Mod(1,p));
}
E = ellrandom(128); ellsea(E, E.a1.mod)
\\VERBOSE flag set -> print diagnostics
E = ellrandom(160); ellsea(E, E.a1.mod, 1)
\\VERBOSE flag set     -> print diagnostics
\\EARLY_ABORT flag set -> abort if a small prime factor is detected
E = ellrandom(192); ellsea(E, E.a1.mod, 1, 1)
\\search a good curve for cryptographic application of size 128bits
ellcrypto(128)
