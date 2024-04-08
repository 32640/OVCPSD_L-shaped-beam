
function [D,Dstr,Tcgm,dTcgmdtheta]= computeElasticityMatrix(C,theta)

n1 = cos(theta);    m1 = sin(theta);
n2 = -sin(theta);   m2 = cos(theta);

Ts = [ n1^2,  m1^2,   -2*m1*n1;
    n2^2,  m2^2,   -2*m2*n2;
    m1*n1, -n1*m1,  n1^2-m1^2];
cta=theta;
nn=cos(cta);mm=sin(cta);

Tcgm=[nn^2, mm^2, 2*nn*mm;
    mm^2, nn^2,-2*nn*mm;
    -nn*mm, nn*mm,nn^2-mm^2];

dTcgmdtheta=[-2*mm*nn,2*mm*nn,2*(nn^2-mm^2);
    2*mm*nn,-2*mm*nn,-2*(nn^2-mm^2);
    -(nn^2-mm^2),nn^2-mm^2,-4*mm*nn];
D = Ts * C * Ts';
Dstr=Tcgm*D;
end