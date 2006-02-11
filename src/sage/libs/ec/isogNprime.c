#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void isogN11()
{GEN jinv=(GEN) TWCURVE[13];

 cv(1,qtwist(avecfromc4c6(stoi(16),stoi(-152)),TWPROD));
 cv(2,qtwist(avecfromc4c6(stoi(496),stoi(20008)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(375376),stoi(229985128)),TWPROD));
 ISOG=25; if (gequal(jinv,gdiv(stoi(-4096),stoi(11)))) PLACE=1;
 else if (gequal(jinv,gdiv(stoi(-122023936),stoi(161051)))) PLACE=2;
 else PLACE=3;
}

void isogN17()
{GEN jinv=(GEN) TWCURVE[13];

 cv(1,qtwist(avecfromc4c6(stoi(33),stoi(-81)),TWPROD));
 cv(2,qtwist(avecfromc4c6(stoi(273),stoi(4455)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(33),stoi(12015)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(4353),stoi(287199)),TWPROD));
 ISOG=4; if (gequal(jinv,gdiv(stoi(35937),stoi(17)))) PLACE=1;
 else if (gequal(jinv,gdiv(stoi(20346417),stoi(289)))) PLACE=2;
 else PLACE=3;
 if (gequal(jinv,gdiv(gmul(stoi(4353),stoi(18948609)),stoi(17))))
 {cv(4,qtwist(avecfromc4c6(stoi(33),stoi(12015)),TWPROD));
  cv(3,qtwist(avecfromc4c6(stoi(4353),stoi(287199)),TWPROD));
 }
}

void isogN19()
{GEN jinv=(GEN) TWCURVE[13];

 cv(1,qtwist(avecfromc4c6(stoi(-32),stoi(8)),TWPROD));
 cv(2,qtwist(avecfromc4c6(stoi(448),stoi(10088)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(36928),stoi(7096328)),TWPROD));
 ISOG=9; if (gequal(jinv,gdiv(stoi(32768),stoi(19)))) PLACE=1;
 else if (gequal(jinv,gdiv(stoi(-89915392),stoi(6859)))) PLACE=2; else PLACE=3;
}

void isogN37()
{GEN jinv=(GEN) TWCURVE[13];

 if (gequal(jinv,gdiv(stoi(110592),stoi(37))))
 {ISOG=1; PLACE=1; cv(1,CURVE); return;}
 cv(1,qtwist(avecfromc4c6(stoi(160),stoi(-2008)),TWPROD));
 cv(2,qtwist(avecfromc4c6(stoi(1120),stoi(36296)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(89920),stoi(26964008)),TWPROD));
 ISOG=9; if (gequal(jinv,gdiv(stoi(4096000),stoi(37)))) PLACE=1;
 else if (gequal(denom(jinv),stoi(37))) PLACE=3; else PLACE=2;
}

void isogNsz()
{GEN jinv=(GEN) TWCURVE[13]; GEN B;

 B=gsub(TWCOND,gmul2n(gun,4));
 if (gequal(jinv,gdiv(gmul(gsqr(B),B),TWCOND)))
 {cv(1,CURVE); cv(2,qtwist(szcN(racine(gsub(TWCOND,gmul2n(gun,6)))),TWPROD));
  ISOG=2; PLACE=1; return;
 }
 B=gsub(gmul2n(gun,8),COND);
 if (gequal(jinv,gdiv(gmul(gsqr(B),B),gsqr(TWCOND))))
 {cv(2,CURVE); cv(1,qtwist(szcP(racine(gsub(TWCOND,gmul2n(gun,6)))),TWPROD));
  ISOG=2; PLACE=2; return;
 }
 ISOG=1; PLACE=1;
}

GEN szcP(GEN u)
{if (gequal(gmod(u,gmul2n(gun,2)),stoi(3))) u=gneg(u);
 return(mineqfromc4c6(gadd(gsqr(u),stoi(48)),
		      gmul(gneg(u),gadd(gsqr(u),stoi(72)))));
}

GEN szcN(GEN u)
{if (gequal(gmod(u,gmul2n(gun,2)),stoi(3))) u=gneg(u);
 return(mineqfromc4c6(gsub(gsqr(u),stoi(192)),
		      gmul(gneg(u),gadd(gsqr(u),stoi(576)))));
}

