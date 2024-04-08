function [StrIJ,StrIJ0,theat,sgema_max,sgema_min]=Stress_T2(DSTR0,xPhys,penal_stress,edofMat,U,nelx,nely,dx1,dy1,act)
nele=nelx*nely;
[x0,y0]=deal(0,0);
[aa,bb]=deal(dx1/2,dy1/2);
sgema=zeros(3,nely*nelx);

for i=1:nele
    Ue=reshape(U(edofMat(i,:),:)',[],1);
    [S]=DB0(x0,y0,aa,bb,cell2mat(DSTR0(i,1)));
    SGM0=xPhys(i).^penal_stress*S*Ue;
    sgema(:,i)=SGM0;
end
CTA=0*sgema(1,:);
sgema_max=zeros(1,nely*nelx);sgema_min=zeros(1,nely*nelx);
sgema_max(act)=(sgema(1,act)+sgema(2,act))/2+(((sgema(1,act)-sgema(2,act))/2).^2+(sgema(3,act)).^2).^0.5;
sgema_min(act)=(sgema(1,act)+sgema(2,act))/2-(((sgema(1,act)-sgema(2,act))/2).^2+(sgema(3,act)).^2).^0.5;
CTA(act)=atan(-2*sgema(3,act)./(sgema(1,act)-sgema(2,act)))./2;
theat = zeros(nele,1);
StrIJ=(1+sign(sgema_max(:)+sgema_min(:)));
[StrIJ]=check(nelx,nely,3,StrIJ,xPhys);
StrIJ=StrIJ(:)';
StrIJ(find(StrIJ>1))=2;
StrIJ(find(StrIJ<=1))=0;
StrIJ0=2-StrIJ;
theat(act)=(CTA(act).*(1+sign(sgema(1,act)-sgema(2,act)))/2+(CTA(act)+pi/2).*(1-sign(sgema(1,act)-sgema(2,act)))/2).*StrIJ(act)/2+...
    (CTA(act).*(1+sign(sgema(1,act)-sgema(2,act)))/2+(CTA(act)+pi/2).*(1-sign(sgema(1,act)-sgema(2,act)))/2+pi/2).*StrIJ0(act)/2;
theat=-theat;


