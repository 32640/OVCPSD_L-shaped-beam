clear; clc
nelx=80; nely=80;nels=40;
dx=1; dy=1;penal=3; penal_stress=0.8;volfrac=0.3; rmin=3;Emin=1e-9;
nele = nelx*nely;pAgg=24;
ndof = 2*(nelx+1)*(nely+1);
C=[2 0.3 0
    0.3 1 0
    0 0 0.25];
%%
n=2;  %% The parameter "n" is represented as "m" in the paper.
x = repmat(0.5,nely,nelx);
theat=repmat(0.5,nely,nelx);
theat_max=zeros(nele,1);
[H,Hs]=filter2d(nelx,nely,rmin);
[iK,jK,edofMat]=iKjK_func(nelx,nely);
%%
iF=[(nely+1)*(nelx-4)-nels:nely+1:(nely+1)*(nelx+1)-nels];
F = sparse(2*iF,1,-1/6,2*(nely+1)*(nelx+1),1);
Lamuda_KS= zeros(2*(nely+1)*(nelx+1),1);
top=[1:(nely+1):(nely+1)*nels+1];
fixeddofs = [2*top 2*top-1];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
U = zeros(ndof,1);
%% 
m=2;
move=0.3;
a0=1;a = zeros(m,1);cmma = 10000*ones(m,1);d= zeros(m,1);
loop = 0; change = 1;
eta=0.5;beta=1;
loopbeta=0;looppAgg=0;looppN=0;
%%
elNrs = reshape(1:nele,nely,nelx);
[pasS] = deal(elNrs(1:nely-nels,nels+1:nelx));
[pasS1] = deal(elNrs(nely-nels+1:nely-nels+2,nelx-4:nelx));
act = setdiff((1:nele)',pasS(:));
x(pasS) = 1e-3; x(pasS1) = 1;
xPhys=x;
while change > 0.001 && loop <500
    loop = loop + 1;
    loopbeta = loopbeta+1;
    looppAgg = looppAgg+1;
    looppN = looppN+1;
    xTilde = x;
    xTilde(:) = (H*x(:))./Hs;xPhysH=0*xTilde(:);
    xPhysH(:) = Heaviside(xTilde,beta,eta);
    xPhys(act)=xPhysH(act);xPhys(pasS1) = 1;
    xval=[x(:);theat(:)];
    xmax=0*xval;
    xmin=0*xval;
    xmax(:) = min(1,xval+move);
    xmin(:) = max(0.001,xval-move);
    if loop==1
        low=0.001*ones(2*nele,1);
        upp=1*ones(2*nele,1);
        xold1=xval;
        xold2=xval;
    end
    xPhys =reshape(xPhys,nely,nelx);
    [Stress0,c,dc,dc_theta,fval_stress,DStressDx,DStressDtheat,theat_max0]=FEA(n,C,xPhys,theat,theat_max,dx,dy,nelx,nely,nele,pAgg,Emin,penal,penal_stress,iK,jK,edofMat,F,freedofs,act);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dxPhys =HeavisideSns(xTilde,beta,eta);
    dv = ones(nely,nelx)/(nele*volfrac);
    dc= H*(dxPhys(:).*dc(:))./Hs;
    dv= H*(dxPhys(:).*dv(:))./Hs;
    DStressDx= H*(dxPhys(:).*DStressDx(:))./Hs;

    [dV,dC,DStressDX,DStressDX2] = deal(zeros(nele,1));
    dC(act)=dc(act);
    dV(act)=dv(act);
    DStressDX(act)=DStressDx(act);

    % 
    v=sum(xPhys(:));M1=1e0;M2=1e1;C_max=400;
    f0val=M1*c/C_max;
    df0dx =M1*[dC(:);dc_theta(:)]/C_max;
    df0dx2=0*df0dx;
    fval=[v/(volfrac*nele)-1;M2*(fval_stress-1)];
    dfdx=zeros(m,2*nele);
    dfdx(1,:)=[dV(:);0*dV(:)];
    dfdx(2,:) =M2*[DStressDX(:);DStressDtheat(:)];
    dfdx2=0*dfdx;

    [xmma,ymma,zmma,lam,xsi,eta0,mu,zet,s,low,upp] = ...
        mmasub(m,2*nele,loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,cmma,d);
    % 
    if loop>2
        xold2 = xold1;
        xold1 = xval;
    elseif loop>1
        xold1 = xval;
    end
    x =xmma(1:nele);
    theat =xmma(nele+1:end);
    x =reshape(x,nely,nelx);
    Obj(loop)=c;
    Con(loop)=mean(x(:));
    theat_max=theat_max0;


    THEAT(loop,:)=theat(:);
    THEAT0(loop,:)=mean(theat(act));
    change = max(abs(xmma(1:nele)-xPhys(:)));

    fprintf('Iter:%3i Obj.:%6.3f fval1.:%6.4f fval2.:%6.4f  change:%6.3f\n',loop,f0val,fval(1),fval(2),change);
    %% 
    figure(2)
    colormap(gray); imagesc(1-xPhys);axis equal; axis tight; axis off; drawnow;
    %%%%%%%%%%%%%%%%
    if  beta < 32  && (loopbeta >= 50 )
        beta = 2*(beta);
        loopbeta = 0;
    end
    if loop>50 && looppAgg >= 50 && pAgg < 100
        pAgg=pAgg+8;
        looppAgg = 0;
    end
end
save('beat1_120_100_0.4_200.mat')
