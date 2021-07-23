clc
clear
close all
[A,B,C,D] = linearize_2wbr(8);
%% linear 2wbr&motor ode45
global i;
global K;

K = [ -51.8355, -54.7619, -17.7365, -6.51892,  0.752264, -0.170483, -587.332,  31.3724;
      -51.8355, -54.7619, -17.7365, -6.51892, -0.170483,  0.752264,  31.3724, -587.332];
init = [0 0 0.01 0 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:1;
i=1;
[t,X] = ode45(@(t,x) linear_ode_2wbr(t,x,A,B),tspan,init,options);
u=-K*X(:,:)';
plot(t,u(1,:));figure
subplot(3,2,1);plot(t,X(:,1)');title('X');
subplot(3,2,2);plot(t,X(:,2));title('dX');
subplot(3,2,3);plot(t,X(:,3)/pi*180);title('TETA');
subplot(3,2,4);plot(t,X(:,4));title('dTETA');
subplot(3,2,5);plot(t,X(:,5)/pi*180);title('w');
subplot(3,2,6);plot(t,X(:,6));title('dw');

function [dx,u] = linear_ode_2wbr(t,x,A,B)
global K;
dx = x;
u = -K*x;
dx = A*x+B*u;
end