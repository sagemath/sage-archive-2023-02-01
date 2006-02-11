#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ellcurv.h"

GEN tatecurve(GEN B,GEN C)
{GEN v=cgetg(6,t_VEC); v[5]=(long) gzero; v[4]=(long) gzero;
 v[3]=(long) gneg(B); v[2]=(long) gneg(B); v[1]=(long) gsub(gun,C);
 return(ellinit0(v,1,ELLACC));
}

GEN cm7(void) {return(mineqfromc4c6(stoi(105),stoi(1323)));}
GEN cm11a(void) {return(mineqfromc4c6(stoi(1441),stoi(54703)));}
GEN cm11b(void) {return(mineqfromc4c6(stoi(352),stoi(-6776)));}
GEN cm11c(void) {return(mineqfromc4c6(stoi(121),stoi(5203)));}
GEN cm19(void) {return(mineqfromc4c6(stoi(1824),stoi(-77976)));}
GEN cm43(void) {return(mineqfromc4c6(stoi(41280),stoi(-8387064)));}
GEN cm67(void) {return(mineqfromc4c6(stoi(353760),stoi(-210408408)));}
GEN cm163(void) {return(mineqfromc4c6(stoi(104372160),
				  gmul(stoi(-40133016),stoi(26569))));}

GEN i2(GEN A)
{GEN C,C4,C6; if (gequal(A,gzero)) return(gzero);
 if (gequalgs(A,-64)) return(gzero);
 C=gaddgs(A,64); C6=gmul(gsubgs(A,8),gsqr(C));
 C4=gmul(gaddgs(A,16),C); return(mineqfromc4c6(C4,C6));
}
GEN i3(GEN A)
{GEN C,C4,C6; if (gequal(A,gzero)) return(gzero);
 if (gequalgs(A,-27)) return(gzero);
 C=gaddgs(A,27); C6=gmul(C,gadd(gaddgs(gmulgs(A,18),-27),gsqr(A)));
 C4=gmul(gaddgs(A,3),C); return(mineqfromc4c6(C4,C6));
}
GEN i4(GEN A)
{if (gequal(A,gmul2n(gun,-1))) return(gzero);
 A=gdiv(gsqr(A),gadd(gmul2n(A,1),gun)); if (gcmp0(A)) return(gzero);
 return(i2(gdivsg(64,A)));
}
GEN i5(GEN A)
{GEN C,C4,C6; if (gequal(A,gzero)) return(gzero);
 C=gadd(gaddgs(gmulgs(A,22),125),gsqr(A));
 C6=gmul(gsqr(C),gadd(gaddgs(gmulgs(A,4),-1),gsqr(A)));
 C4=gmul(C,gadd(gaddgs(gmulgs(A,10),5),gsqr(A)));
 return(mineqfromc4c6(C4,C6));
}
GEN i6(GEN A)
{if (gequalgs(A,-4)) return(gzero);
 A=gdiv(gmul(gsqr(A),gaddgs(gmul2n(A,1),9)),gmulgs(gsqr(gaddgs(A,4)),27));
 if (gcmp0(A)) return(gzero); return(i3(gdivsg(27,A)));
}
GEN i7(GEN A)
{GEN C,C4,C6; if (gequal(A,gzero)) return(gzero);
 C=gadd(gaddgs(gmulgs(A,13),49),gsqr(A));
 C6=gmul(gsqr(A),gadd(gaddgs(gmulgs(A,14),63),gsqr(A)));
 C6=gmul(C,gadd(C6,gaddgs(gmulgs(A,70),-7)));
 C4=gmul(C,gadd(gaddgs(gmulgs(A,5),1),gsqr(A)));
 return(mineqfromc4c6(C4,C6));
}
GEN i8(GEN A)
{if (gequalgs(A,-4)) return(gzero); if (gcmp0(A)) return(gzero);
 A=gdiv(gdeux,gmul(A,gaddgs(A,4))); return(i4(A));
}
GEN i9(GEN A)
{if (gcmp0(A)) return(gzero);
 A=ginv(gmul(A,gadd(gsqr(A),gmulsg(3,gadd(A,gun)))));
 if (gcmp0(A)) return(gzero); return(i3(gdivsg(27,A)));
}
GEN i10(GEN A)
{GEN C=gaddgs(A,2);
 if (gcmp0(C)) return(gzero); if (gcmp0(A)) return(gzero);
 A=gdiv(gaddgs(gmul2n(A,1),5),gmul(A,gmul(C,gsqr(gsqr(C)))));
 if (gcmp0(A)) return(gzero); return(i2(gdivsg(64,A)));
}
GEN i12(GEN A)
{if (gequalgs(A,-6)) return(gzero); if (gcmp0(A)) return(gzero);
 A=gdivsg(36,gmul(A,gaddgs(A,6))); return(i6(A));
}
GEN i13(GEN A)
{GEN C1,C2,C4,C6; if (gcmp0(A)) return(gzero);
 C1=gadd(gsqr(A),gaddgs(gmulgs(A,5),13));
 C2=gadd(gsqr(A),gaddgs(gmulgs(A,6),13));
 C4=gmul(gsqr(A),gadd(gsqr(A),gaddgs(gmulgs(A,7),20)));
 C4=gmul(gmul(C1,C2),gadd(C4,gaddgs(gmulgs(A,19),1)));
 C6=gmul(gsqr(gsqr(A)),gadd(gsqr(A),gaddgs(gmulgs(A,10),46)));
 C6=gadd(C6,gmul(gsqr(A),gaddgs(gmulgs(A,108),122)));
 C6=gmul(gmul(C1,gsqr(C2)),gadd(C6,gaddgs(gmulgs(A,38),-1)));
 return(mineqfromc4c6(C4,C6));
}
GEN i15a() {return(mineqfromc4c6(stoi(25),stoi(1475)));}
GEN i15b() {return(mineqfromc4c6(stoi(145),stoi(-2105)));}
GEN i16(GEN A) {return(i8(gmul2n(gmul(A,gaddgs(A,4)),-1)));}
GEN i17a() {return(mineqfromc4c6(stoi(145945),stoi(-55755325)));}
GEN i17b() {return(mineqfromc4c6(stoi(31705),stoi(6328675)));}
GEN i18(GEN A)
{if (gcmp0(A)) return(gzero);
 A=gdivsg(36,gmul(A,gadd(gsqr(A),gaddgs(gmulgs(A,6),12)))); return(i6(A));
}
GEN i21a() {return(mineqfromc4c6(stoi(225),stoi(-3537)));}
GEN i21b() {return(mineqfromc4c6(stoi(-135),stoi(243)));}
GEN i25(GEN A)
{GEN P,Q,R,C4,C6,C; if (gcmp0(A)) return(gzero);
 C=gadd(gsqr(A),gaddgs(gmulgs(A,2),5));
 Q=gmul(gsqr(A),gadd(gsqr(A),gaddgs(gmulgs(A,5),15)));
 Q=gadd(Q,gaddgs(gmulgs(A,25),25));
 P=gmul(gsqr(gsqr(gsqr(A))),gadd(gsqr(A),gaddgs(gmulgs(A,10),55)));
 P=gadd(P,gmul(gmul(gsqr(gsqr(A)),gsqr(A)),gaddgs(gmulgs(A,200),525)));
 P=gadd(P,gmul(gsqr(gsqr(A)),gaddgs(gmulgs(A,1010),1425)));
 P=gadd(P,gmul(gsqr(A),gaddgs(gmulgs(A,1400),875)));
 P=gadd(P,gaddgs(gmulgs(A,250),5)); C4=gmul(P,C);
 R=gsub(P,gmulgs(gadd(gmul(A,Q),gun),6)); R=gmul(gsqr(C),R);
 C6=gmul(gsqr(A),gadd(gsqr(A),gaddgs(gmulgs(A,4),9)));
 C6=gadd(C6,gaddgs(gmulgs(A,10),5)); C6=gmul(C6,R);
 return(mineqfromc4c6(C4,C6));
}

GEN i27() {return(mineqfromc4c6(gzero,stoi(-216)));}
GEN i37() {return(mineqfromc4c6(stoi(385),stoi(-8225)));}

GEN sc(GEN A)
{GEN C,D;
 C=cgetg(6,t_VEC); C[1]=(long) gzero; C[2]=(long) gcopy(A);
 C[3]=(long) gzero; C[4]=(long) gneg(gadd(A,stoi(3))); C[5]=(long) gun;
 D=ellinit0(C,1,ELLACC); return(mineqfromc4c6((GEN) D[10],(GEN) D[11]));
}

GEN specialcurve(char *INP)
{char which[16]; GEN P,Q,T,U,C4,C6,CF,A,B,V,E,C; int i,n,OK;

 if (INP[0]=='j') {INP++; return(jcurve(INP));}
 if (!strcmp(INP,"sn") || !strcmp(INP,"SN"))
 {scanf("%s",which);
  if (!strcmp(which,"2a") || !strcmp(which,"2A"))
  {scanf("%s",which); U=lltoi(atoll(which));
   C4=gadd(gsqr(U),stoi(48)); C6=gmul(U,gadd(gsqr(U),stoi(72)));
   if (atoll(which)==15) printf("SN2A CURVE\n");
   if (atoll(which)==0) printf("SN2A CURVE\n");
   if (gequal(gmodulcp(U,gsqr(gdeux)),stoi(3)))
   {if (isprime(gadd(gmul2n(gun,6),gsqr(U)))) printf("SN2A CURVE\n");}
   return(mineqfromc4c6(C4,C6));
  }
  if (!strcmp(which,"2b") || !strcmp(which,"2B"))
  {scanf("%s",which); U=lltoi(atoll(which));
   C4=gadd(gsqr(gmul2n(U,2)),stoi(48));
   C6=gmul(gmul2n(U,2),gadd(gsqr(gmul2n(U,2)),stoi(72)));
   if (isprime(gadd(gmul2n(gun,2),gsqr(U)))) printf("SN2B CURVE\n");
   if (atoll(which)==11) printf("SN2B CURVE\n");
   if (atoll(which)==-11) printf("SN2B CURVE\n");
   if (atoll(which)==-2) printf("SN2B CURVE\n");
   if (atoll(which)==2) printf("SN2B CURVE\n");
   if (atoll(which)==0) printf("SN2B CURVE\n");
   return(mineqfromc4c6(C4,C6));
  }
  if (!strcmp(which,"2c") || !strcmp(which,"2C"))
  {scanf("%s",which); U=lltoi(atoll(which));
   C4=gadd(gsqr(U),gmul2n(gadd(U,gun),4));
   C6=gmul(gadd(U,gmul2n(gun,3)),
	   gadd(gsqr(U),gsub(gmul2n(U,4),gmul2n(gun,3))));
   if (gequal(gmodulcp(U,gsqr(gdeux)),gun)) return(mineqfromc4c6(C4,C6));
   P=isprimepower(gabs(U,-1)); Q=isprimepower(gabs(gadd(U,gmul2n(gun,4)),-1));
   if (gequal(U,gneg(gun))) printf("SN2C CURVE\n");
   if (gequal(U,stoi(-17))) printf("SN2C CURVE\n");
   if (!gequal(P,gzero) && !gequal(Q,gzero))
   {if ((gequal(gmodulcp(P,gsqr(gdeux)),stoi(3))) ||
	(gequal(gmodulcp(Q,gsqr(gdeux)),stoi(3)))) printf("SN2C CURVE\n");
   }
   return(mineqfromc4c6(C4,C6));
  }
  if (!strcmp(which,"3a") || !strcmp(which,"3A"))
  {scanf("%s",which); U=lltoi(atoll(which));
   C4=gmul(gadd(U,stoi(3)),gadd(gadd(gmul(U,gsqr(U)),gmul(gsqr(U),stoi(9))),
				gadd(gmul(stoi(27),U),stoi(3))));
   C6=gadd(gmul(U,gadd(gmul(U,gadd(U,stoi(18))),stoi(135))),stoi(504));
   C6=gadd(gmul(U,gadd(gmul(U,gadd(gmul(C6,U),stoi(891))),stoi(486))),
	   stoi(-27)); C6=gneg(C6);
   if (!gequal(isprimepower(gadd(gsqr(U),
				gadd(gmul(U,stoi(9)),stoi(27)))),gzero))
   {CF=factor(U); n=itos((GEN) matsize(CF)[1]); OK=1;
    for(i=1;i<=n;i++)
    {T=(GEN) ((GEN) CF[1])[i]; if (gequal(gmodulcp(T,stoi(3)),gun)) OK=0;}
    if (OK) printf("SN3A CURVE\n");
   }
   return(mineqfromc4c6(C4,C6));
  }
 }
 if (!strcmp(INP,"tw") || !strcmp(INP,"TW"))
 {scanf("%s",which); U=lltoi(atoll(which));
  scanf("%s",INP); return(qtwist(getcurve(INP),U));
 }

 if (!strcmp(INP,"t2") || !strcmp(INP,"T2"))
 {scanf("%s",which); A=gp_read_str(which); scanf("%s",which); B=gp_read_str(which);
  if (gequal(B,gzero)) return(gzero);
  if (gequal(gsqr(A),gmul2n(B,2))) return(gzero);
  V=cgetg(6,t_VEC); V[5]=(long) gzero; V[4]=(long) B;
  V[3]=(long) gzero; V[2]=(long) A; V[1]=(long) gzero;
  E=ellinit0(V,1,ELLACC); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t3") || !strcmp(INP,"T3"))
 {scanf("%s",which); A=gp_read_str(which); scanf("%s",which); B=gp_read_str(which);
  if (gequal(B,gzero)) return(gzero);
  if (gequal(gmul(A,gsqr(A)),gmul(B,stoi(27)))) return(gzero);
  V=cgetg(6,t_VEC); V[5]=(long) gzero; V[4]=(long) gzero;
  V[3]=(long) B; V[2]=(long) gzero; V[1]=(long) A;
  E=ellinit0(V,1,ELLACC); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t4") || !strcmp(INP,"T4"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(gneg(A),gmul2n(gun,-4))) return(gzero);
  E=tatecurve(A,gzero); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t22") || !strcmp(INP,"T22"))
 {scanf("%s",which); A=gp_read_str(which); scanf("%s",which); B=gp_read_str(which);
  V=cgetg(6,t_VEC); V[5]=(long) gzero; V[4]=(long) gmul(A,B);
  V[3]=(long) gzero; V[2]=(long) gadd(A,B); V[1]=(long) gzero;
  E=ellinit0(V,1,ELLACC); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t5") || !strcmp(INP,"T5"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  E=tatecurve(A,A); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t6") || !strcmp(INP,"T6"))
 {scanf("%s",which); A=gp_read_str(which);
  if (gequal(A,gzero) || gequal(A,gneg(gun)) || gequal(A,gdiv(gun,stoi(-9))))
    return(gzero); E=tatecurve(gadd(A,gsqr(A)),A);
  return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t7") || !strcmp(INP,"T7"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(A,gun)) return(gzero); B=gsub(gsqr(A),A);
  E=tatecurve(gmul(A,B),B); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t8") || !strcmp(INP,"T8"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(A,gun) || gequal(A,gmul2n(gun,-1))) return(gzero);
  B=gmul(gsub(A,gun),gsub(gmul2n(A,1),gun));
  E=tatecurve(B,gdiv(B,A)); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t24") || !strcmp(INP,"T24"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(gabs(A,-1),gmul2n(gun,-2))) return(gzero);
  E=tatecurve(gsub(gsqr(A),gmul2n(gun,-4)),gzero);
  return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t9") || !strcmp(INP,"T9"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(A,gun)) return(gzero);
  B=gadd(gun,gmul(A,gsub(A,gun))); A=gmul(A,gsub(B,gun));
  E=tatecurve(gmul(A,B),A); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t10") || !strcmp(INP,"T10"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(A,gun) || gequal(A,gmul2n(gun,-1))) return(gzero);
  B=gdiv(gsqr(A),gsub(A,gsqr(gsub(A,gun)))); A=gmul(A,gsub(B,gun));
  E=tatecurve(gmul(A,B),A); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t12") || !strcmp(INP,"T12"))
 {scanf("%s",which); A=gp_read_str(which); if (gequal(A,gzero)) return(gzero);
  if (gequal(A,gun) || gequal(A,gmul2n(gun,-1))) return(gzero);
  B=gdiv(gsub(gmul(stoi(3),gmul(A,gsub(gun,A))),gun),gsub(A,gun));
  V=gadd(B,A); A=gmul(gdiv(B,gsub(gun,A)),gsub(V,gun));
  E=tatecurve(gmul(A,V),A); return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t26") || !strcmp(INP,"T26"))
 {scanf("%s",which); A=gp_read_str(which);
  if (gequal(gsub(gabs(A,-1),gun),gdeux)) return(gzero);
  if (gequal(gsub(A,gun),gsqr(gdeux))) return(gzero);
  if (gequal(gsub(A,gun),gmul2n(gun,3))) return(gzero);
  B=gdiv(gsub(stoi(10),gmul2n(A,1)),gsub(gsqr(A),stoi(9)));
  E=tatecurve(gmul(B,gadd(B,gun)),B);
  return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"t28") || !strcmp(INP,"T28"))
 {scanf("%s",which); A=gp_read_str(which);
  if (gequal(A,gzero)) return(gzero);
  if (gequal(gneg(A),gmul2n(gun,-1))) return(gzero);
  if (gequal(gneg(A),gmul2n(gun,-2))) return(gzero);
  B=gdiv(gmul(A,gadd(gdeux,gmul2n(A,3))),gsub(gmul2n(gsqr(A),3),gun));
  A=gmul(gsub(gmul2n(B,1),gun),gsub(B,gun)); E=tatecurve(A,gdiv(A,B));
  return(mineqfromc4c6((GEN) E[10],(GEN) E[11]));
 }
 if (!strcmp(INP,"i2") || !strcmp(INP,"I2"))
 {scanf("%s",which); A=gp_read_str(which); return(i2(A));}
 if (!strcmp(INP,"i3") || !strcmp(INP,"I3"))
 {scanf("%s",which); A=gp_read_str(which); return(i3(A));}
 if (!strcmp(INP,"i4") || !strcmp(INP,"I4"))
 {scanf("%s",which); A=gp_read_str(which); return(i4(A));}
 if (!strcmp(INP,"i5") || !strcmp(INP,"I5"))
 {scanf("%s",which); A=gp_read_str(which); return(i5(A));}
 if (!strcmp(INP,"i6") || !strcmp(INP,"I6"))
 {scanf("%s",which); A=gp_read_str(which); return(i6(A));}
 if (!strcmp(INP,"i7") || !strcmp(INP,"I7"))
 {scanf("%s",which); A=gp_read_str(which); return(i7(A));}
 if (!strcmp(INP,"i8") || !strcmp(INP,"I8"))
 {scanf("%s",which); A=gp_read_str(which); return(i8(A));}
 if (!strcmp(INP,"i9") || !strcmp(INP,"I9"))
 {scanf("%s",which); A=gp_read_str(which); return(i9(A));}
 if (!strcmp(INP,"i10") || !strcmp(INP,"I10"))
 {scanf("%s",which); A=gp_read_str(which); return(i10(A));}
 if (!strcmp(INP,"i11") || !strcmp(INP,"I11")) return(cm11a());
 if (!strcmp(INP,"i11a") || !strcmp(INP,"I11A")) return(cm11a());
 if (!strcmp(INP,"i11b") || !strcmp(INP,"I11B")) return(cm11b());
 if (!strcmp(INP,"i11c") || !strcmp(INP,"I11C")) return(cm11c());
 if (!strcmp(INP,"i12") || !strcmp(INP,"I12"))
 {scanf("%s",which); A=gp_read_str(which); return(i12(A));}
 if (!strcmp(INP,"i13") || !strcmp(INP,"I13"))
 {scanf("%s",which); A=gp_read_str(which); return(i13(A));}
 if (!strcmp(INP,"i14") || !strcmp(INP,"I14")) return(cm7());
 if (!strcmp(INP,"i15") || !strcmp(INP,"I15")) return(i15a());
 if (!strcmp(INP,"i15a") || !strcmp(INP,"I15A")) return(i15a());
 if (!strcmp(INP,"i15b") || !strcmp(INP,"I15B")) return(i15b());
 if (!strcmp(INP,"i16") || !strcmp(INP,"I16"))
 {scanf("%s",which); A=gp_read_str(which); return(i16(A));}
 if (!strcmp(INP,"i17") || !strcmp(INP,"I17")) return(i17a());
 if (!strcmp(INP,"i17a") || !strcmp(INP,"I17A")) return(i17a());
 if (!strcmp(INP,"i17b") || !strcmp(INP,"I17B")) return(i17b());
 if (!strcmp(INP,"i18") || !strcmp(INP,"I18"))
 {scanf("%s",which); A=gp_read_str(which); return(i18(A));}
 if (!strcmp(INP,"i19") || !strcmp(INP,"I19")) return(cm19());
 if (!strcmp(INP,"i21") || !strcmp(INP,"I21")) return(i21a());
 if (!strcmp(INP,"i21a") || !strcmp(INP,"I21A")) return(i21a());
 if (!strcmp(INP,"i21b") || !strcmp(INP,"I21B")) return(i21b());
 if (!strcmp(INP,"i25") || !strcmp(INP,"I25"))
 {scanf("%s",which); A=gp_read_str(which); return(i25(A));}
 if (!strcmp(INP,"i27") || !strcmp(INP,"I27")) return(i27());
 if (!strcmp(INP,"i37") || !strcmp(INP,"I37")) return(i37());
 if (!strcmp(INP,"i43") || !strcmp(INP,"I43")) return(cm43());
 if (!strcmp(INP,"i67") || !strcmp(INP,"I67")) return(cm67());
 if (!strcmp(INP,"i163") || !strcmp(INP,"I163")) return(cm163());
 if (!strcmp(INP,"sc") || !strcmp(INP,"SC"))
 {scanf("%s",which); A=gp_read_str(which); return(sc(A));}
 return(gzero);
}


