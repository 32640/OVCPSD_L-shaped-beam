load beat1_120_100_0.4_1.mat
theat=pi*(theat-0.5)./n;
theat0=theat(:)+theat_max(:);
cta=-theat0(:);
xxx=xPhys;
TTT=reshape( xxx, nely, nelx );
TTT(find(TTT>(0.5)))=1;
TTT(find(TTT<0.5))=0;
%%%%%%%%%%%%%%%%%
sgema100 =reshape(Stress0,nely,nelx);
sgema1000 =TTT.*sgema100 ;
[A1,A2]=max(sgema1000(:))
[B1,B2]=min(sgema1000(:))

sgema1000(find(sgema1000==0))=500;

figure(1)
colormap(jet); imagesc(sgema1000); axis equal; axis off; drawnow;


[x,y]=meshgrid(0:1:nelx,0:1:nely);
[x1,y1]=meshgrid(0.5:1:nelx,0.5:1:nely);

ll=0.4;
x2=x1(:)+ll*cos(cta);
y2=y1(:)+ll*sin(cta);
x1=x1(:)-ll*cos(cta);
y1=y1(:)-ll*sin(cta);


x1=TTT(:).*x1(:);
y1=TTT(:).*y1(:);
x2=TTT(:).*x2(:);
y2=TTT(:).*y2(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(x1)
    line([x1(i)+0.5,x2(i)+0.5],[y1(i)+0.5,y2(i)+0.5],'color','w','linewidth',1)
end



