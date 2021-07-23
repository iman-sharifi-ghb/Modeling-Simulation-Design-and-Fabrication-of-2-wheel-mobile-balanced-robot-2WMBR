%clc 
clear 
close all 
%% nonlinear 2wbr&motor ode45
init = [0 0 0.00001 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
[t,X] = ode45(@(t,x) ode_2wbr(t,x),tspan,init,options);
subplot(2,2,1);plot(t,X(:,1));title('X');
subplot(2,2,2);plot(t,X(:,2));title('Xdot');
subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('TETA');
subplot(2,2,4);plot(t,X(:,4));title('TETAdot');

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

% Equations
dx = x;
TR = 0;
TL = 0;

% assumed equations-------------------------------------------------
n1 = Mb+2*Mw+2*J/r^2;
n2 = Mb * l*cos(x(3));
n3 = -Mb*l*(x(6)^2+x(4)^2)*sin(x(3))+2/r*Ca*(x(2)/r-x(4))-(TR+TL)/r;
n4 = Mb*l*cos(x(3));
n5 = (I2+Mb*l^2);
n6 = (I3 - I1 - Mb*l^2)*x(6)^2*sin(x(3))*cos(x(3))-Mb*l*g*sin(x(3))-2*Ca*(x(2)/r-x(4))+(TL+TR);
n7 = (I3 + 2*K + Mw*d^2/2+J*d^2/r^2/2) -(I3 - I1 - Mb*l^2)*(sin(x(3)))^2;
n8 = (Mb*l*x(2)-2*(I3 - I1 - Mb*l^2)*x(4)*cos(x(3)))*x(6)*sin(x(3))+Ca*x(6)*d^2/r^2/2-(TR-TL)*d/2/r;

F1 = (n2*n6-n3*n5)/(n1*n5-n2*n4);
F2 = (n3*n4-n1*n6)/(n1*n5-n2*n4);
F3 = -n8/n7;

dx = [x(2); F1; x(4); F2; x(6); F3];
end
