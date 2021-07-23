clc
clear
close all
[A,B,C,D] = linearize_2wbr(4);
%%
sys  = ss(A,B,C,D);
sys = tf(sys)
% rlocus(sys);figure
% nyquist(sys);figure
% bode(sys);figure
p0 = vpa(pole(sys),8)
z0 = vpa(zero(sys),8)
co = ctrb(A,B);
rank(co,1)
% state feedback
p =[ -5
 -154.19679
  -10.6733926
 -20.7081816]
K = place(A,B,p);
disp('K=')
disp(vpa(K,6))
Acl = A-B*K;
eig(Acl);
sys1  = ss(Acl,B,C,D);
step(1/3.14*180*sys1);figure
sys1 = tf(sys1);
% rlocus(sys1);figure
% nyquist(sys1);figure
% bode(sys1)
% LQR 
Q = diag([1 1 100 1]);
R = 0.1;
N = [0 0 0 0]';
[KK,S,e] = lqr(A,B,Q,R,N);
AA = A-B*KK;
vpa(eig(AA),6)
sys2 = ss(AA,B,C,D);
step(1/3.14*180*sys2);%figure
sys2 = tf(sys2);
% rlocus(sys2);figure
% nyquist(sys2);figure
% bode(sys2)





















