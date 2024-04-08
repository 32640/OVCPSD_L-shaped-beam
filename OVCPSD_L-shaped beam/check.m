%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,dc,xPhys)
 dc =reshape(dc,nely,nelx);
dcn=zeros(nely,nelx);xPhys(:,:)=1;
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*xPhys(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(xPhys(j,i)*sum);
  end
end