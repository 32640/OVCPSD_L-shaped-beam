function [SK0,SK1,Kee,DSTR,DSTR0,DStrDcta00]=lk_H8(DH,theat,dx,dy,nele,Emin,xPhys,penal,n)
theat=theat(:);
SK0=cell(nele,1);SK1=cell(nele,1);Kee=cell(nele,1);
DSTR=cell(nele,1);DSTR0=cell(nele,1);
DStrDcta00=cell(nele,1);
for e=1:nele
 [D,Dstr,Tcgm,dTcgmdtheta]= computeElasticityMatrix(DH,theat(e));
 [dDdtheat,dStrDdtheta]= computedDdtheta1(DH,D,theat(e),n,Tcgm,dTcgmdtheta);
 DSTR{e,1}= Dstr;DSTR0{e,1}= D;
 DStrDcta00{e,1}= dStrDdtheta;
 [K,dKdtheat]=Stiffness(dx,dy,D,dDdtheat);      
 Kee{e,1}= K;
 sK=(Emin+xPhys(e).^penal*(1-Emin)).*(K(:));
 SK0{e,1}=sK;
 SK1{e,1}=dKdtheat;
end




