clc
clear
close all
[A,B,C,D] = linearize_2wbr(4);
%% linear 2wbr&motor ode45
u =4.753;
init = [0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:0.5;
[t,X] = ode45(@(t,x) linear_ode_2wbr(t,x,A,B,u),tspan,init,options);
subplot(2,2,1);plot(t,X(:,1));title('X');
subplot(2,2,2);plot(t,X(:,2));title('dX');
subplot(2,2,3);plot(t,X(:,3)/pi*180);title('TETA');
subplot(2,2,4);plot(t,X(:,4));title('dTETA');

function dx = linear_ode_2wbr(t,x,A,B,u)
dx = x;
dx = A*x+B*u;
end