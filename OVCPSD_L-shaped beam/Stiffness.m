function [K,dKdtheat]=Stiffness(dx,dy,D,dDdtheat)
node=[0 0;dx 0;dx dy;0 dy];ts=0.004;
point_1=1/sqrt(3);
gauss_point(1,1)=-point_1;  gauss_point(1,2)=-point_1; 
gauss_point(2,1)=point_1;   gauss_point(2,2)=-point_1; 
gauss_point(3,1)=point_1;   gauss_point(3,2)=point_1;
gauss_point(4,1)=-point_1;  gauss_point(4,2)=point_1;
K=0; dKdtheat=0;
for i=1:4
    s=gauss_point(i,1); t=gauss_point(i,2);
    [J,B]=Matrix(s,t,node);
    K=K+B'*D*B*det(J)*ts;
    dKdtheat=dKdtheat+B'*dDdtheat*B*det(J)*ts;
end