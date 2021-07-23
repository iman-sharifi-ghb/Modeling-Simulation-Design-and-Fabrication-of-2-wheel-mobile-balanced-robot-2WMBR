%clc 
clear 
close all 
%% nonlinear 2wbr&motor ode45
init = [0 0 0 0 0 0 0 0 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
[t,X] = ode45(@(t,x) ode_2wbr(t,x),tspan,init,options);
subplot(2,2,1);plot(t,X(:,1));title('X');
subplot(2,2,2);plot(t,X(:,2));title('Xdot');
subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('TETA');
subplot(2,2,4);plot(t,X(:,4));title('TETAdot');figure
plot(t,X(:,3)/3.14*180);title('TETA')

function dx = ode_2wbr(t,x)
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
Km = p.Km;
Kb = p.Kb;
R = p.R;
L = p.L;
J_e = 2e-3;
B = p.J_e;
% assumed equations-------------------------------------------------
n1 = Mb+2*Mw+2*J/r^2;
n2 = I2 + Mb*l^2;
n3 = I3 - I1 - Mb*l^2;
n4 = I3 + 2*K + Mw*d^2/2;
% Equations
u1 = 4.753;
u2 = 4.753;
dx = x;
TR = Km*x(11);
TL = Km*x(12);
F1 = 1/(n1-1/n2*(Mb*l*cos(x(3)))^2)*((TR+TL)/r+...
    (Mb*l*(x(6)^2+x(4)^2)*sin(x(3))-2/r*Ca*(x(2)/r-x(4)))-...
    Mb/n2*l*cos(x(3))*(-n3*x(6)^2*sin(x(3))*cos(x(3))+...
    Mb*l*g*sin(x(3))+2*Ca*(x(2)/r-x(4))))+...
    Mb*l/n2*cos(x(3))*(TR+TL);
F2 = 1/n2*(-(TR+TL)+(-n3*x(6)^2*sin(x(3))*cos(x(3))+...
    Mb*l*g*sin(x(3))+2*Ca*(x(2)/r-x(4)))-...
    Mb*l*cos(x(3))*(1/(n1-1/n2*(Mb*l*cos(x(3)))^2)*((TR+TL)/r+...
    (Mb*l*(x(6)^2+x(4)^2)*sin(x(3))-2/r*Ca*(x(2)/r-x(4)))-...
    Mb/n2*l*cos(x(3))*(-n3*x(6)^2*sin(x(3))*cos(x(3))+...
    Mb*l*g*sin(x(3))+2*Ca*(x(2)/r-x(4))))+Mb*l/n2*cos(x(3))*(TR+TL)));
F3 = 1/(n4+J*d^2/2/r^2-n3*(sin(x(3))^2))*((TL-TR)*d/2/r-...
    (Mb*l*x(2)-2*n3*x(4)*cos(x(3)))*x(6)*sin(x(3))-Ca*x(6)*d^2/(2*r^2));

dx = [x(2); F1; x(4); F2; x(6); F3; x(9); x(10);
    (Km*x(11)-B*x(9)+Ca*x(4))/J_e;
    (Km*x(12)-B*x(10)+Ca*x(4))/J_e;
    (-Kb/L)*x(9)-R/L*x(11)+1/L*u1;
    (-Kb/L)*x(10)-R/L*x(12)+1/L*u2];
end
