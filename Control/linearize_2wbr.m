function [A,B,C,D] = linearize_2wbr(u)
u = int8(u);
% 2wbr parameters---------------------------------------------------
p = IIwbr_parameters();
d = p.d;
l = p.l;
r = p.r;
Mb = p.Mb;
Mw = p.Mw;
g = p.g;
I1 = p.I1;
I2 = p.I2;
I3 = p.I3;
J = p.J;
K = p.K;
Ca = p.Ca;
% dc motor parameters--------------------------------------------------
Km = 1.34;
Kb = 1;
R = 1.7;
L = 0.05;
% J_m = 2.403e-5;
% J1 = 1.5842e-5;
% N = 59;
% J_e = J1+N^2*J_m;
J_e = 2e-3;
B = p.B;
% assumed equations-------------------------------------------------
n1 = Mb+2*Mw+2*J/r^2;
n2 = I2 + Mb*l^2;
n3 = I3 - I1 - Mb*l^2;
n4 = I3 + 2*K + Mw*d^2/2;

if u == 12
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 u1 u2
    % x1=x, x2=dx x3=teta, x4=dteta, x5=sai, x6=dsai, x7=tetaR, x8=tetaL,
    % x9=wR, x10=wL, x11=iR, x12=iL, u1=vR, u2=vL
    TR = Km*x11;
    TL = Km*x12;
    F1 = 1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+...
        Mb*l/n2*cos(x3)*(TR+TL);
    F2 = 1/n2*(-(TR+TL)+(-n3*x6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4))-...
        Mb*l*cos(x3)*(1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+Mb*l/n2*cos(x3)*(TR+TL)));
    F3 = 1/(n4+J*d^2/2/r^2-n3*(sin(x3)^2))*((TL-TR)*d/2/r-...
        (Mb*l*x2-2*n3*x4*cos(x3))*x6*sin(x3)-Ca*x6*d^2/(2*r^2));
    
    F = [x2; F1; x4; F2; x6; F3; x9; x10;
        (Km*x11-B*x9+Ca*x4)/J_e;
        (Km*x12-B*x10+Ca*x4)/J_e;
        (-Kb/L)*x9-R/L*x11+1/L*u1;
        (-Kb/L)*x10-R/L*x12+1/L*u2];
    F = simplify(F);
    
    x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];
    u = [u1;u2];
    
    A = jacobian(F,x);
    A = simplify(A);
    B = jacobian(F,u);
    B = simplify(B);
    
    A = subs(A,[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,u1,u2],...
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,u1,u2],...
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [0 0 1 0 0 0 0 0];
    D = zeros(1,2);
elseif u==8
    syms x1 x2 x3 x4 x5 x6 x7 x8 u1 u2
    % x1=x, x2=dx x3=teta, x4=dteta, x5=wR, x6=wL, x7=iR, x8=iL,u1=vR,u2=vL
    x_6 = 0;
    TR = Km*x7;
    TL = Km*x8;
    F1 = 1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+...
        Mb*l/n2*cos(x3)*(TR+TL);
    F2 = 1/n2*(-(TR+TL)+(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4))-...
        Mb*l*cos(x3)*(1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+Mb*l/n2*cos(x3)*(TR+TL)));
    
    F = [x2; F1; x4; F2;
        (TR-B*x5+Ca*x4)/J_e;
        (TL-B*x6-Ca*x4)/J_e;
        (-Kb/L)*x5-R/L*x7+1/L*u1;
        (-Kb/L)*x6-R/L*x8+1/L*u2];
    F = simplify(F);
    
    x = [x1;x2;x3;x4;x5;x6;x7;x8];
    u = [u1;u2];
    % state space matrixs
    A = jacobian(F,x);
    A = simplify(A);
    B = jacobian(F,u);
    B = simplify(B);
    
    A = subs(A,[x1,x2,x3,x4,x5,x6,x7,x8,u1,u2],...
        [0,0,0,0,0,0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,x5,x6,x7,x8,u1,u2],...
        [0,0,0,0,0,0,0,0,0,0]);
    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [0 0 1 0 0 0 0 0];
    D = zeros(1,2);
elseif u==6
    syms x1 x2 x3 x4 x5 x6 u
    % x1=x, x2=dx x3=teta, x4=dteta, x5=w, x6=i,u=v
    x_6 = 0;
    T = Km*x6;
    TR = T;
    TL = T;
    F1 = 1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+...
        Mb*l/n2*cos(x3)*(TR+TL);
    F2 = 1/n2*(-(TR+TL)+(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4))-...
        Mb*l*cos(x3)*(1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+Mb*l/n2*cos(x3)*(TR+TL)));
    
    F = [x2; F1; x4; F2;
        (T-B*x5+Ca*x4)/J_e;
        (-Kb/L)*x5-R/L*x6+1/L*u];
    F = simplify(F);
    
    x = [x1;x2;x3;x4;x5;x6];
    % state space matrixs
    A = jacobian(F,x);
    A = simplify(A);
    B = jacobian(F,u);
    B = simplify(B);
    
    A = subs(A,[x1,x2,x3,x4,x5,x6,u],...
        [0,0,0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,x5,x6,u],...
        [0,0,0,0,0,0,0]);
    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [0 0 1 0 0 0];
    D = zeros(1,1);
elseif u==4
    syms x1 x2 x3 x4 u
    % x1=teta, x2=dteta, x3=w, x4=i,u=v
    x_6 = 0;% SAI is not important

    F1 = 1/(n1-1/n2*(Mb*l*cos(x3))^2)*((2*u)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+...
        Mb*l/n2*cos(x3)*(2*u);
    F2 = 1/n2*(-(2*u)+(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4))-...
        Mb*l*cos(x3)*(1/(n1-1/n2*(Mb*l*cos(x3))^2)*((2*u)/r+...
        (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
        Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
        Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+Mb*l/n2*cos(x3)*(2*u)));
    
    F = [x2;F1;
         x4;F2];
    F = simplify(F);
    
    x = [x1;x2;x3;x4];
    % state space matrixs
    A = jacobian(F,x);
    A = simplify(A);
    B = jacobian(F,u);
    B = simplify(B);
    
    A = subs(A,[x1,x2,x3,x4,u],...
        [0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,u],...
        [0,0,0,0,0]);
    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [0 0 1 0];
    D = zeros(1,1);
end
%     else
%         % stepper motor parameters
%         Km = 1;
%         Nr = 25;
%         R = 3;
%         L = 0.1;
%         Rm = 10;
%         B = 0.003;
%         J = 1e-6;
%         syms x1 x2 x3 x4 x5 x6 x7 x8 u1 u2
%         % x1=x, x2=dx x3=teta, x4=dteta, x5=iA, x6=iB, x7=w,
%         % x8=teta,u1=vA,u2=vB
%         x_6 = 0;
%         TR = Km*x7;
%         TL = Km*x8;
%         F1 = 1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
%             (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
%             Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
%             Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+...
%             Mb*l/n2*cos(x3)*(TR+TL);
%         F2 = 1/n2*(-(TR+TL)+(-n3*x_6^2*sin(x3)*cos(x3)+...
%             Mb*l*g*sin(x3)+2*Ca*(x2/r-x4))-...
%             Mb*l*cos(x3)*(1/(n1-1/n2*(Mb*l*cos(x3))^2)*((TR+TL)/r+...
%             (Mb*l*(x_6^2+x4^2)*sin(x3)-2/r*Ca*(x2/r-x4))-...
%             Mb/n2*l*cos(x3)*(-n3*x_6^2*sin(x3)*cos(x3)+...
%             Mb*l*g*sin(x3)+2*Ca*(x2/r-x4)))+Mb*l/n2*cos(x3)*(TR+TL)));
%
%         eA = -Km*x7*sin(Nr*x8);
%         eB = Km*x7*sin(Nr*x8);
%         Td = Ca*(x7-x4);
%         Te = -Km*(x5-eA/Rm)*sin(Nr*x8)+Km*(x6-eB/Rm)*cos(Nr*x8)-Td*sin(4*Nr*x8);
%
%         F = [x2; F1; x4; F2;
%             (u1-R*x5-eA)/L;
%             (u2-R*x6-eB)/L;
%             1/J*(Te-B*x7);
%             x7];
%         F = simplify(F);
%
%         x = [x1;x2;x3;x4;x5;x6;x7;x8];
%         u = [u1;u2];
%         % state space matrixs
%         A = jacobian(F,x);
%         A = simplify(A);
%         B = jacobian(F,u);
%         B = simplify(B);
%
%         A = subs(A,[x1,x2,x3,x4,x5,x6,x7,x8,u1,u2],...
%             [0,0,0,0,0,0,0,0,0,0]);
%         B = subs(B,[x1,x2,x3,x4,x5,x6,x7,x8,u1,u2],...
%             [0,0,0,0,0,0,0,0,0,0]);
%         A = vpa(A,6);
%         B = vpa(B,6);
%         A = double(A);
%         B = double(B);
%         C = [0 0 1 0 0 0 0 0];
%         D = zeros(1,2);
%     end
end