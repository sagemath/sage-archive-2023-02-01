#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

void checkexotic()
{GEN jinv=(GEN) CURVE[13]; GEN edisc=(GEN) TWCURVE[12];

 if (gequal(jinv,gzero)) CM=-3;
 if (gequal(jinv,stoi(1728))) CM=-1; if (gequal(jinv,stoi(8000))) CM=-2;

 if (gequal(jinv,gzero) && (gequal(edisc,stoi(-27))))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=3;
   exo27(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=1; exo27();}
 }
 if (gequal(jinv,stoi(-12288000)) && (gequal(edisc,stoi(-243))))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=4;
   exo27(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=2; exo27();}
 }

 if (gequal(jinv,gmul2n(stoi(3375),-1)))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=2;
   exo21A(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=1; exo21B();}
 }
 if (gequal(jinv,gmul2n(stoi(-140625),-3)))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=2;
   exo21B(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=1; exo21A();}
 }
 if (gequal(jinv,gmul2n(gmul(stoi(-5745),stoi(33005025)),-7)))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=4;
   exo21A(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=3; exo21B();}
 }
 if (gequal(jinv,gmul2n(stoi(-1159088625),-21)))
 {if (gequal(gmod(TWPROD,stoi(3)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-3)); PLACE=4;
   exo21B(); TWPROD=gmul(TWPROD,stoi(-3));
  }
  else {PLACE=3; exo21A();}
 }

 if (gequal(jinv,stoi(-3375)))
 {if (gequal(gmod(TWPROD,stoi(7)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-7)); PLACE=3;
   exo14(); TWPROD=gmul(TWPROD,stoi(-7));
  }
  else {PLACE=1; exo14();}
 }
 if (gequal(jinv,stoi(16581375)))
 {if (gequal(gmod(TWPROD,stoi(7)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-7)); PLACE=4;
   exo14(); TWPROD=gmul(TWPROD,stoi(-7));
  }
  else {PLACE=2; exo14();}
 }

 if (gequal(jinv,gmul2n(stoi(-25),-1)))
 {if (gequal(gmod(TWPROD,stoi(5)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(5)); PLACE=3;
   exo15B(); TWPROD=gmul(TWPROD,stoi(5));
  }
  else {PLACE=1; exo15A();}
 }
 if (gequal(jinv,gmul2n(stoi(-349938025),-3)))
 {if (gequal(gmod(TWPROD,stoi(5)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(5)); PLACE=4;
   exo15B(); TWPROD=gmul(TWPROD,stoi(5));
  }
  else {PLACE=2; exo15A();}
 }
 if (gequal(jinv,gmul2n(stoi(-121945),-5)))
 {if (gequal(gmod(TWPROD,stoi(5)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(5)); PLACE=3;
   exo15A(); TWPROD=gmul(TWPROD,stoi(5));
  }
  else {PLACE=1; exo15B();}
 }
 if (gequal(jinv,gmul2n(stoi(46969655),-15)))
 {if (gequal(gmod(TWPROD,stoi(5)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(5)); PLACE=4;
   exo15A(); TWPROD=gmul(TWPROD,stoi(5));
  }
  else {PLACE=2; exo15B();}
 }

 if (gequal(jinv,stoi(-24729001)))
 {if (gequal(gmod(TWPROD,stoi(11)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-11)); PLACE=2;
   exo11C(); TWPROD=gmul(TWPROD,stoi(-11));
  }
  else {PLACE=1; exo11A();}
 }
 if (gequal(jinv,stoi(-32768)))
 {if (gequal(gmod(TWPROD,stoi(11)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-11)); PLACE=2;
   exo11B(); TWPROD=gmul(TWPROD,stoi(-11));
  }
  else {PLACE=1; exo11B();}
 }
 if (gequal(jinv,stoi(-121)))
 {if (gequal(gmod(TWPROD,stoi(11)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-11)); PLACE=2;
   exo11A(); TWPROD=gmul(TWPROD,stoi(-11));
  }
  else {PLACE=1; exo11C();}
 }

 if (gequal(jinv,gmul2n(stoi(-297756989),-1)))
 {if (gequal(gmod(TWPROD,stoi(17)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(17)); PLACE=2;
   exo17B(); TWPROD=gmul(TWPROD,stoi(17));
  }
  else {PLACE=1; exo17A();}
 }
 if (gequal(jinv,gmul2n(stoi(-882216989),-17)))
 {if (gequal(gmod(TWPROD,stoi(17)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(17)); PLACE=2;
   exo17A(); TWPROD=gmul(TWPROD,stoi(17));
  }
  else {PLACE=1; exo17B();}
 }

 if (gequal(jinv,stoi(-9317))) {exo37(); PLACE=1;}
 if (gequal(jinv,gmul(stoi(-1997597),gsqr(stoi(285371))))) {exo37(); PLACE=2;}

 if (gequal(jinv,stoi(-884736)))
 {if (gequal(gmod(TWPROD,stoi(19)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-19)); PLACE=2;
   exo19(); TWPROD=gmul(TWPROD,stoi(-19));
  }
  else {PLACE=1; exo19();}
 }
 if (gequal(jinv,stoi(-884736000)))
 {if (gequal(gmod(TWPROD,stoi(43)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-43)); PLACE=2;
   exo43(); TWPROD=gmul(TWPROD,stoi(-43));
  }
  else {PLACE=1; exo43();}
 }
 if (gequal(jinv,gmul(stoi(-5280),gsqr(stoi(5280)))))
 {if (gequal(gmod(TWPROD,stoi(67)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-67)); PLACE=2;
   exo67(); TWPROD=gmul(TWPROD,stoi(-67));
  }
  else {PLACE=1; exo67();}
 }
 if (gequal(jinv,gmul(stoi(-640320),gsqr(stoi(640320)))))
 {if (gequal(gmod(TWPROD,stoi(163)),gzero))
  {TWPROD=gdiv(TWPROD,stoi(-163)); PLACE=2;
   exo163(); TWPROD=gmul(TWPROD,stoi(-163));
  }
  else {PLACE=1; exo163();}
 }
}

void exo11A()
{cv(1,qtwist(avecfromc4c6(stoi(1441),stoi(54703)),TWPROD)); ISOG=11; CM=-11;
 cv(2,qtwist(avecfromc4c6(stoi(14641),stoi(-6925193)),TWPROD));
}

void exo11B()
{cv(1,qtwist(avecfromc4c6(stoi(352),stoi(-6776)),TWPROD)); ISOG=11; CM=-11;
 cv(2,qtwist(avecfromc4c6(stoi(42592),stoi(9018856)),TWPROD));
}

void exo11C()
{cv(1,qtwist(avecfromc4c6(stoi(121),stoi(5203)),TWPROD)); ISOG=11; CM=-11;
 cv(2,qtwist(avecfromc4c6(stoi(174361),stoi(-72809693)),TWPROD));
}

void exo14()
{cv(1,qtwist(avecfromc4c6(stoi(105),stoi(1323)),TWPROD)); ISOG=14; CM=-7;
 cv(2,qtwist(avecfromc4c6(stoi(1785),stoi(75411)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(5145),stoi(-453789)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(87465),stoi(-25865973)),TWPROD));
}

void exo15A()
{cv(1,qtwist(avecfromc4c6(stoi(25),stoi(1475)),TWPROD)); ISOG=15;
 cv(2,qtwist(avecfromc4c6(stoi(6025),stoi(467675)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(3625),stoi(-263125)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(-26375),stoi(1941875)),TWPROD));
}

void exo15B()
{cv(1,qtwist(avecfromc4c6(stoi(145),stoi(-2105)),TWPROD)); ISOG=15;
 cv(2,qtwist(avecfromc4c6(stoi(-1055),stoi(15535)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(625),stoi(184375)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(150625),stoi(58459375)),TWPROD));
}

void exo17A()
{cv(1,qtwist(avecfromc4c6(stoi(145945),stoi(-55755325)),TWPROD)); ISOG=17;
 cv(2,qtwist(avecfromc4c6(stoi(9162745),
			  gmul(stoi(2088025),stoi(14891))),TWPROD));
}

void exo17B()
{cv(1,qtwist(avecfromc4c6(stoi(31705),stoi(6328675)),TWPROD)); ISOG=17;
 cv(2,qtwist(avecfromc4c6(stoi(42178105),
			  gmul(stoi(-35496425),stoi(7717))),TWPROD));
}

void exo19()
{cv(1,qtwist(avecfromc4c6(stoi(1824),stoi(-77976)),TWPROD)); ISOG=19; CM=-19;
 cv(2,qtwist(avecfromc4c6(stoi(658464),stoi(534837384)),TWPROD));
}

void exo21A()
{cv(1,qtwist(avecfromc4c6(stoi(225),stoi(-3537)),TWPROD)); ISOG=21;
 cv(2,qtwist(avecfromc4c6(stoi(-1215),stoi(-6561)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(4545),stoi(622431)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(465345),stoi(317440863)),TWPROD));
}

void exo21B()
{cv(1,qtwist(avecfromc4c6(stoi(-135),stoi(243)),TWPROD)); ISOG=21;
 cv(2,qtwist(avecfromc4c6(stoi(2025),stoi(95499)),TWPROD));
 cv(3,qtwist(avecfromc4c6(stoi(51705),stoi(-11757069)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(40905),stoi(-16805637)),TWPROD));
}

void exo27()
{cv(1,qtwist(avecfromc4c6(gzero,stoi(-216)),TWPROD)); ISOG=27;
 cv(2,qtwist(avecfromc4c6(stoi(1440),stoi(-54648)),TWPROD));
 cv(3,qtwist(avecfromc4c6(gzero,stoi(5832)),TWPROD));
 cv(4,qtwist(avecfromc4c6(stoi(12960),stoi(1475496)),TWPROD));
}

void exo37()
{cv(1,qtwist(avecfromc4c6(stoi(385),stoi(-8225)),TWPROD)); ISOG=37;
 cv(2,qtwist(avecfromc4c6(stoi(9987985),
			  gmul(stoi(2758525),stoi(11443))),TWPROD));
}

void exo43()
{GEN A=qtwist(avecfromc4c6(stoi(41280),stoi(-8387064)),TWPROD);
 cv(1,A); cv(2,qtwist(A,stoi(-43))); ISOG=43; CM=-43;
}

void exo67()
{GEN A=qtwist(avecfromc4c6(stoi(353760),stoi(-210408408)),TWPROD);
 cv(1,A); cv(2,qtwist(A,stoi(-67))); ISOG=67; CM=-67;
}

void exo163()
{GEN A=qtwist(avecfromc4c6(stoi(104372160),
			   gmul(stoi(-40133016),stoi(26569))),TWPROD);
 cv(1,A); cv(2,qtwist(A,stoi(-163))); ISOG=163; CM=-163;
}
