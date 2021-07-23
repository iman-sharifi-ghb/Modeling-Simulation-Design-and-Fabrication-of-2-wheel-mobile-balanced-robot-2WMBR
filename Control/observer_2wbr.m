clc
clear
close all
[A,B,C,D] = linearize_2wbr(4);
%% control 
sys  = ss(A,B,C,D);% sys = tf(sys);
p0 = vpa(pole(sys),8)
z0 = vpa(zero(sys),8)
co = ctrb(A,B);
Controlability = size(A,1)-rank(co,1);
if Controlability == 0
    fprintf('system is controllable\n');
else
    fprintf('system is not controllable\n');
end
ob = obsv(A,C);
Observability = ctrb(A',C') - (obsv(A,C))';
if Observability == 0
    fprintf('system is observable');
else
    fprintf('system is not observable');
end

ControllerPoles=[-1
  -270.0555
 -6.7006074
   -6.944848];
K=place(A,B,ControllerPoles);

ObserverPoles = [-1
  -3
 -2
   -10];
L=place(A',C',ObserverPoles).';

Aaug=[A    -B*K
      L*CC   A-B*K-L*CC];
  
Baug=[B
      B];

Caug=[CC -DD*K];

Daug=DD;

H=ss(Aaug,Baug,Caug,Daug);
H = tf(H);
step(H(3,1))

% H=H/dcgain(H);
% 
% tf=25;
% t=linspace(0,tf,1000)';
% 
% r=1*(t>=0);%-1*(t>=5)+3*(t>=10)-2*(t>=15);
% 
% x0=[0 0 0 0 0 0 -0.5 0 zeros(1,8)];
% 
% [y, ~, xaug]=lsim(H,r,t,x0);
% disp(xaug)

% x=xaug(:,[1 2]);
% xhat=xaug(:,[3 4]);
% ex=x-xhat;
% 
% figure;
% step(sys,H);

% figure;
% subplot(2,1,1);
% plot(t,x(:,1),'b',t,xhat(:,1),'r');
% legend('x_1','Estimation of x_1');
% xlabel('t');
% ylabel('x_1');
% subplot(2,1,2);
% plot(t,ex(:,1),'k');
% legend('Estimation Error of x_1');
% xlabel('t');
% ylabel('e_{x_1}');
% 
% figure;
% subplot(2,1,1);
% plot(t,x(:,2),'b',t,xhat(:,2),'r');
% legend('x_2','Estimation of x_2');
% xlabel('t');
% ylabel('x_2');
% subplot(2,1,2);
% plot(t,ex(:,2),'k');
% legend('Estimation Error of x_2');
% xlabel('t');
% ylabel('e_{x_2}');

