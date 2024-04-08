function [Stress0,c,dc,dc_theta,fval_stress,DStressDx,DStressDtheat,theat_max0]=FEA(n,DH,xPhys,theat,theat_max,dx,dy,nelx,nely,nele,pAgg,Emin,penal,penal_stress,iK,jK,edofMat,F,freedofs,act)
theat=pi*(theat-0.5)./n;
theat0=theat(:)+theat_max(:);
U=0*F;
[SK0,SK1,Kee,DSTR,DSTR0,DStrDcta00]=lk_H8(DH,theat0,dx,dy,nele,Emin,xPhys,penal,n);
K = sparse(iK(:),jK(:),cell2mat(SK0)); K = (K+K')/2;
U(freedofs)=K(freedofs,freedofs)\F(freedofs);
%%
[StrIJ,StrIJ0,theat_max0,sgema_max,sgema_min]=Stress_T2(DSTR0,xPhys,penal_stress,edofMat,U,nelx,nely,dx,dy,act);
[Stress0,fval_stress,gama,H0,Htheat,CP]=Stress_T1(StrIJ,StrIJ0,DSTR,xPhys,penal_stress,edofMat,U,nelx,nely,pAgg,dx,dy,DStrDcta00,act,sgema_max,sgema_min);
[c,dc,dc_theta,DStressDx,DStressDtheat]=Obj_Sen(CP,gama,xPhys,edofMat,U,K,nelx,nely,freedofs,penal,Kee,SK1,Emin,H0,Htheat,act);
