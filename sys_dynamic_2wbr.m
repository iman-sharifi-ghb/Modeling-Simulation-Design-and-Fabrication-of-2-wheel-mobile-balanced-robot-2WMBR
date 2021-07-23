function dZ=AM_dynamic(t,Z)

dZ=Z;

global Mb Mw I1 I2 I3 K J l d r g

q=Z(1:6); dq=Z(7:12);
q1=q(1);q2=q(2);q3=q(3); q4=q(4);q5=q(5);q6=q(6);
dq1=dq(1); dq2=dq(2); dq3=dq(3); dq4=dq(4);dq5=dq(5);dq6=dq(6);

M = Mfunc(I1,I2,I3,J,K,Mb,Mw,d,q3,r);
B = Bfunc(I1,I3,Mb,dq4,g,l,q3);
a = afunc(d,q4,r);
adot = adotfunc(dq4,q4);

F=-B;

M=[M -a';-a zeros(3)];

F=[F;adot*dq];

X=M\F;

dZ=[dq;X(1:6)];
landa(1:3) = X(7:9);
% disp(landa)
