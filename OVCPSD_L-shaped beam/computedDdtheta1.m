function [dDdtheta,dStrDdtheta]= computedDdtheta1(DH,D,theta,n,Tcgm,dTcgmdtheta)

n1 = cos(theta);    m1 = sin(theta);
n2 = -sin(theta);   m2 = cos(theta);

Ts = [ n1^2,  m1^2,   -2*m1*n1;
    n2^2,  m2^2,   -2*m2*n2;
    m1*n1, -n1*m1,  n1^2-m1^2];


dTsdtheta=[-2*cos(theta)*sin(theta),     2*cos(theta)*sin(theta), 2*sin(theta)^2 - 2*cos(theta)^2
    2*cos(theta)*sin(theta),    -2*cos(theta)*sin(theta), 2*cos(theta)^2 - 2*sin(theta)^2
    cos(theta)^2 - sin(theta)^2, sin(theta)^2 - cos(theta)^2,        -4*cos(theta)*sin(theta)];

dTstdtheta=[-2*cos(conj(theta))*sin(conj(theta)),         2*cos(conj(theta))*sin(conj(theta)), cos(conj(theta))^2 - sin(conj(theta))^2
    2*cos(conj(theta))*sin(conj(theta)),        -2*cos(conj(theta))*sin(conj(theta)), sin(conj(theta))^2 - cos(conj(theta))^2
    2*sin(conj(theta))^2 - 2*cos(conj(theta))^2, 2*cos(conj(theta))^2 - 2*sin(conj(theta))^2,    -4*cos(conj(theta))*sin(conj(theta))];


dDdtheta = pi*(dTsdtheta*DH*Ts'+Ts*DH*dTstdtheta)./n;
dStrDdtheta = pi*dTcgmdtheta*D./n+Tcgm*dDdtheta;

end