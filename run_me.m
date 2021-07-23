clc
close all
clear
%%
addpath('Generate');

global Mb Mw I1 I2 I3 K J l d r g Ca

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

dt=0.01;
timeSpan=0:dt:10;
tsize=length(timeSpan);

q_0=[0;0;0;0;0;0];
dq_0=[0;0;0;0;0;0];
init=[q_0;dq_0;0];
mu_0=[0;0;0];

options=odeset('maxstep',1e-3);
%%
t0=clock;
[~,zIMM]=ode45(@IMM_dynamic,timeSpan,[q_0;dq_0;mu_0],options);
t1=clock;
timeIMM=etime(t1,t0);
disp(['IMM sim time: ' num2str(timeIMM) '(s)'])
%%
t0=clock;
[t,zAM]=ode45(@AM_dynamic,timeSpan,init,options);
t1=clock;
timeAM=etime(t1,t0);
disp(['AM sim time: ' num2str(timeAM) '(s)'])
%%
subplot(3,2,1);plot(t,zIMM(:,1));title('X');
subplot(3,2,2);plot(t,zIMM(:,7));title('Xdot');
subplot(3,2,3);plot(t,zIMM(:,3));title('TETA');
subplot(3,2,4);plot(t,zIMM(:,9));title('TETAdot');
subplot(3,2,5);plot(t,zIMM(:,4));title('SAI');
subplot(3,2,6);plot(t,zIMM(:,10));title('SAIdot');
figure
subplot(3,2,1);plot(t,zAM(:,13));title('X');
subplot(3,2,2);plot(t,zAM(:,7));title('Xdot');
subplot(3,2,3);plot(t,zAM(:,3));title('TETA');
subplot(3,2,4);plot(t,zAM(:,9));title('TETAdot');
subplot(3,2,5);plot(t,zAM(:,4));title('SAI');
subplot(3,2,6);plot(t,zAM(:,10));title('SAIdot');
%%
IM.q1=zIMM(:,1);   IM.dq1=zIMM(:,7);
IM.q2=zIMM(:,2);   IM.dq2=zIMM(:,8);
IM.q3=zIMM(:,3);   IM.dq3=zIMM(:,9);
IM.q4=zIMM(:,4);   IM.dq4=zIMM(:,10);
IM.q5=zIMM(:,5);   IM.dq5=zIMM(:,11);
IM.q6=zIMM(:,6);   IM.dq6=zIMM(:,12);

AM.q1=zAM(:,1);   AM.dq1=zAM(:,7);
AM.q2=zAM(:,2);   AM.dq2=zAM(:,8);
AM.q3=zAM(:,3);   AM.dq3=zAM(:,9);
AM.q4=zAM(:,4);   AM.dq4=zAM(:,10);
AM.q5=zAM(:,5);   AM.dq5=zAM(:,11);
AM.q6=zAM(:,6);   AM.dq6=zAM(:,12);
%%
%Constraints
CE_IM= CEfunc(d,IM.dq1,IM.dq2,IM.dq4,IM.dq5,IM.dq6,IM.q4,r);
CE_AM= CEfunc(d,AM.dq1,AM.dq2,AM.dq4,AM.dq5,AM.dq6,AM.q4,r);

%Mechanical Energy
E_IM=Efunc(I1,I2,I3,J,K,Mb,Mw,d,IM.dq1,IM.dq2,IM.dq3,IM.dq4,g,l,IM.q3,IM.q4,r);
E_AM=Efunc(I1,I2,I3,J,K,Mb,Mw,d,AM.dq1,AM.dq2,AM.dq3,AM.dq4,g,l,AM.q3,AM.q4,r);
 
ER_IM=(E_IM-E_IM(1))/E_IM(1);
ER_AM=(E_IM-E_AM(1))/E_AM(1);

%% plot
figure
hold on; grid on
plot(timeSpan,ER_IM,'r-')
plot(timeSpan,ER_AM,'g-')
legend('IM','AM')
xlabel('time(s)')
title('Mechanical Energy Error')
figure
hold on; grid on
plot(timeSpan,CE_IM,'r-')
plot(timeSpan,CE_AM,'g-')
legend('IM','AM')
xlabel('time(s)')
title('Constraint Error')
figure
hold on; grid on
plot(timeSpan,IM.q1,'r-','linewidth',4)
plot(timeSpan,AM.q1,'g-','linewidth',3)
legend('IM','AM')
xlabel('time(s)')
title('X')
figure
% hold on; grid on
% plot(timeSpan,IM.q2,'r-','linewidth',4)
% plot(timeSpan,AM.q2,'g-','linewidth',3)
% legend('IM','AM')
% xlabel('time(s)')
% title('q_2')
% figure
hold on; grid on
plot(timeSpan,IM.q3*1/3.14*180,'r-','linewidth',4)
plot(timeSpan,AM.q3*1/3.14*180,'g-','linewidth',3)
legend('IM','AM')
xlabel('time(s)')
title('TETA') 
% figure
% hold on; grid on
% plot(timeSpan,IM.q4,'r-','linewidth',4)
% plot(timeSpan,AM.q4,'g-','linewidth',3)
% legend('IM','AM')
% xlabel('time(s)')
% title('q_4')
%%
function dZ=AM_dynamic(t,Z)

dZ=Z;

global Mb Mw I1 I2 I3 K J l d r g Ca

q=Z(1:6); dq=Z(7:12);
q1=q(1);q2=q(2);q3=q(3); q4=q(4);q5=q(5);q6=q(6);
dq1=dq(1); dq2=dq(2); dq3=dq(3); dq4=dq(4);dq5=dq(5);dq6=dq(6);

M = Mfunc(I1,I2,I3,J,K,Mb,Mw,d,l,q3,q4,r);
B = Bfunc(I1,I3,J,Mb,Mw,dq1,dq2,dq4,g,l,q3,q4,r);
a = afunc(d,q4,r);
adot = adotfunc(dq4,q4);
TL = 0.05;
TR = 0.05;
eq1 = TL-Ca*(dq5-dq3);
eq2 = TR-Ca*(dq6-dq3);
Q = [0;0;-(eq1+eq2);0;eq1;eq2];

F=-B+Q;

M=[M -a';-a zeros(3)];

F=[F;adot*dq];

X=M\F;

dZ=[dq;X(1:6);dq1*cos(q4)+dq2*sin(q4)];
landa(1:3) = double(X(7:9));
% disp(landa)
end
function dZ=IMM_dynamic(t,Z)

dZ=Z;

global Mb Mw I1 I2 I3 K J l d r g Ca

q=Z(1:6); dq=Z(7:12);
q1=q(1);q2=q(2);q3=q(3); q4=q(4);q5=q(5);q6=q(6);
dq1=dq(1); dq2=dq(2); dq3=dq(3); dq4=dq(4);dq5=dq(5);dq6=dq(6);

M = Mfunc(I1,I2,I3,J,K,Mb,Mw,d,l,q3,q4,r);
B = Bfunc(I1,I3,J,Mb,Mw,dq1,dq2,dq4,g,l,q3,q4,r);
a = afunc(d,q4,r);
adot = adotfunc(dq4,q4);
TL = 0.05;
TR = 0.05;
eq1 = TL-Ca*(dq5-dq3);
eq2 = TR-Ca*(dq6-dq3);
Q = [0;0;-(eq1+eq2);0;eq1;eq2];
F=-B+Q;

M_IMM=[eye(6) zeros(6) zeros(6,3);
    zeros(6) M -a';
    zeros(3,6) -a zeros(3)];

F_IMM=[dq;F;adot*dq];


dZ=M_IMM\F_IMM;
end

