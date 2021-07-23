clc 
clear 
close all

syms t
syms d l r 
syms Mb Mw g
syms I1 I2 I3 J K
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6
% q1=xc, q2=yc q3=teta, x4=sai, x5=gamaR, x6=gamaL
 
% vL = [dq1-d/2*dq4*cos(q4);dq2-d/2*dq4*sin(q4);0];
% vR = [dq1+d/2*dq4*cos(q4);dq2+d/2*dq4*sin(q4);0];
% vB = [dq1+l*dq3*cos(q3)*cos(q4)-l*dq4*sin(q3)*sin(q4);
%       dq2+l*dq3*cos(q3)*sin(q4)+l*dq4*sin(q3)*cos(q4);l*dq4*sin(q3)];

C2N = [cos(q4) -sin(q4) 0;
       sin(q4) cos(q4) 0;
       0 0 1];
B2C = [cos(q3) 0 sin(q3);
       0 1 0;
       -sin(q3) 0 cos(q3)];

dx = dq1*cos(q4)+dq2*sin(q4);

vL = [dx - d/2*dq4;0;0];
vR = [dx + d/2*dq4;0;0];
vB = [dx;0;0] + B2C*[0;l*dq3;l*dq4*sin(q3)];

wL = [0;(1/r)*(dx - d/2*dq4);dq4];
wR = [0;(1/r)*(dx + d/2*dq4);dq4];
wB = [-dq4*sin(q3);dq3;dq4*cos(q3)];

v1=[K J K];
v2=[I1 I2 I3];
I_R = diag(v1);
I_L = I_R;
I_B = diag(v2);

T_trans = 1/2*Mw*(vL.'*vL) + 1/2*Mw*(vR.'*vR) + 1/2*Mb*(vB.'*vB);
T_rot = 1/2*(wL.'*I_L*wL) + 1/2*(wR.'*I_R*wR) + 1/2*(wB.'*I_B*wB);
T = T_trans + T_rot;
V = Mb*g*l*cos(q3);
L = T-V;

% constraints
c1 = dx-d/2*dq4-r*dq5;
c2 = dx+d/2*dq4-r*dq6;
c3 = dq1*sin(q4)-dq2*cos(q4);

% Generalized coordinates
q=[q1;q2;q3;q4;q5;q6];
dq=[dq1;dq2;dq3;dq4;dq5;dq6];

%%
dL_dq=jacobian(L,q);
dL_ddq=jacobian(L,dq);

M=jacobian(dL_ddq,dq);
B=diff(dL_ddq.',t)-dL_dq.';
C=[c1;c2;c3];
CE=simplify(C.'*C);

a=jacobian(C,dq);
adot=a;
for n=1:size(a,1)
    adot(n,:)=(jacobian(a(n,:),q)*dq).';   
end

M=simplify(M);
B=simplify(B);
E=simplify(T+V);
a=simplify(a);
adot=simplify(adot);
%%
MFunc=matlabFunction(M,'File','generate\Mfunc');
BFunc=matlabFunction(B,'File','generate\Bfunc');
aFunc=matlabFunction(a,'File','generate\afunc');
adotFunc=matlabFunction(adot,'File','generate\adotfunc');
EFunc=matlabFunction(E,'File','generate\Efunc');
CEFunc=matlabFunction(CE,'File','generate\CEfunc');

