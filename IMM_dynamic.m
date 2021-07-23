function dZ=IMM_dynamic(t,Z)

dZ=Z;

global Mb Mw I1 I2 I3 K J l d r g Ca

q=Z(1:6); dq=Z(7:12);
q1=q(1);q2=q(2);q3=q(3); q4=q(4);q5=q(5);q6=q(6);
dq1=dq(1); dq2=dq(2); dq3=dq(3); dq4=dq(4);dq5=dq(5);dq6=dq(6);

M = Mfunc(I1,I2,I3,J,K,Mb,Mw,d,l,q3,q4,r);
B = Bfunc(I1,I3,J,Mb,dq1,dq2,dq3,dq4,g,l,q3,q4,r);
a = afunc(d,q4,r);
adot = adotfunc(dq4,q4);
TL = 1;
TR = 1;
eq1 = TL-Ca*(dq5-dq3);
eq2 = TR-Ca*(dq6-dq3);
Q = [0;0;-(eq1+eq2);0;eq1;eq2];
F=-B+Q;

M_IMM=[eye(6) zeros(6) zeros(6,3);
    zeros(6) M -a';
    zeros(3,6) -a zeros(3)];

F_IMM=[dq;F;adot*dq];


dZ=M_IMM\F_IMM;

