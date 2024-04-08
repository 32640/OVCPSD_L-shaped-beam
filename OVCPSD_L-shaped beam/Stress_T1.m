function [sgema1,SMS_1_KS,gama,H0,Htheat,CP]=Stress_T1(StrIJ,StrIJ0,DSTR,xPhys,penal_stress,edofMat,U,nelx,nely,p,dx1,dy1,DStrDcta00,act,sgema_max,sgema_min)
nele=nelx*nely;
[x0,y0]=deal(0,0);
[aa,bb]=deal(dx1/2,dy1/2);
sgema=zeros(3,nely*nelx);

for i=1:nele
    Ue=reshape(U(edofMat(i,:),:)',[],1);
    [S]=DB0(x0,y0,aa,bb,cell2mat(DSTR(i,1)));
    SGM0=xPhys(i).^penal_stress*S*Ue;
    sgema(:,i)=SGM0;
end

sgema1=0*sgema(1,:);
sgema1(:)=sgema(1,:);
sgema100 =reshape(sgema1,nely,nelx);
figure(1)
colormap(jet); imagesc(sgema100); axis equal; axis off; drawnow;
SMG_Max_T=120;SMG_Max_C=-100;
SMS_1=0*sgema100(:);dHoff=0*sgema100(:);
SMS_1(act)=sgema1(act)./(2*SMG_Max_T).*StrIJ(act)+sgema1(act)./(2*SMG_Max_C).*StrIJ0(act);
dHoff(act)=1/(2*SMG_Max_T).*StrIJ(act)+1/(2*SMG_Max_C).*StrIJ0(act);
sgm_max=max(sgema1(:))
sgm_min=min(sgema1(:))

%% KS
sd2=0*sgema100;sd3=0*sgema100;dKS_dHoff=0*sgema100;
sd2(act)=exp(p.*(SMS_1(act)));
sd3(act)=sum(exp(p.*(SMS_1(act))));
sd4=log(sum(exp(p.*(SMS_1(act)))));
SMS_1_KS0=sd4/p      
dKS_dHoff(act)=sd2(act)./(sd3(act));
[SMS_1_KS,N_max]=max(SMS_1)
CP=SMS_1_KS/SMS_1_KS0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0=zeros(nely*nelx,3);Htheat=zeros(nely*nelx,3);
gama=zeros(2*(nely+1)*(nelx+1),3); index_matrix=edofMat';
for i=act'
    u=reshape(U(edofMat(i,:),:)',[],1);
    [S0]=DB0(x0,y0,aa,bb,cell2mat(DSTR(i,1)));
    [dSdtheat]=DB0(x0,y0,aa,bb,cell2mat(DStrDcta00(i)));
    H0(i,:)=penal_stress*xPhys(i).^(penal_stress-1)*dKS_dHoff(i)*dHoff(i)*S0*u;
    Htheat(i,:)=xPhys(i).^penal_stress*dKS_dHoff(i)*dHoff(i)*dSdtheat*u;
    index=index_matrix(:,i);
    gama(index,:)=gama(index)+xPhys(i).^penal_stress*dKS_dHoff(i)*dHoff(i)*S0';
end


