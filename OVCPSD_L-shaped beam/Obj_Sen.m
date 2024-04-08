function [c,dc,dc_theta,DStressDx,DStressDtheat]=Obj_Sen(CP,gama,xPhys,edofMat,U,K,nelx,nely,freedofs,penal,Kee,SK1,Emin,H0,Htheat,act)
  nele=nelx*nely;
  Lamuda_KS= zeros(2*(nely+1)*(nelx+1),3);
  Lamuda_KS(freedofs,:)= K(freedofs,freedofs)\gama(freedofs,:);
  Lamuda_KS1=Lamuda_KS(:,1);
  Lamuda_KS1=Lamuda_KS1(edofMat);
  Lamuda_KS2=Lamuda_KS(:,2);
  Lamuda_KS2=Lamuda_KS2(edofMat);
  Lamuda_KS3=Lamuda_KS(:,3);
  Lamuda_KS3=Lamuda_KS3(edofMat);
  Lamuda_KS0=[Lamuda_KS1,Lamuda_KS2,Lamuda_KS3];

    UU=U(edofMat);
    c=0;dc=0*ones(nele,1);dc_theta=0*ones(nele,1);
    Lamuda_x=zeros(nely*nelx,3);Lamuda_theat=zeros(nely*nelx,3);
    for i=act'
        Ue=UU(i,:)';
        Lamuda_KS00=Lamuda_KS0(i,:);
        Lamuda_KSe=reshape(Lamuda_KS00,8,3);
        c=c+xPhys(i)^penal*Ue'*cell2mat(Kee(i))*Ue;
        dc(i,1)=-penal*xPhys(i)^(penal-1)*Ue'*cell2mat(Kee(i))*Ue;
        dc_theta(i,1)=-xPhys(i)^(penal)*Ue'*cell2mat(SK1(i))*Ue;
        Lamuda_x(i,:) =penal*(1-Emin)*xPhys(i).^(penal-1).*Lamuda_KSe'*cell2mat(Kee(i))*Ue;
        Lamuda_theat(i,:) =(1-Emin)*xPhys(i).^(penal).*Lamuda_KSe'*cell2mat(SK1(i))*Ue; 
    end  
    dStressdx=H0-Lamuda_x;
    DStressDx=CP*dStressdx(:,1);
    dStressdtheat=Htheat-Lamuda_theat;
    DStressDtheat=CP*dStressdtheat(:,1);