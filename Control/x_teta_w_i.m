clc
clear
close all
[A,B,C,D] = linearize_2wbr(4);
%%
sys  = ss(A,B,C,D);
sys = tf(sys)
% rlocus(sys);%figure
% % nyquist(sys);figure
% % bode(sys);figure
% sys = c2d(sys,0.0001)
% [zd,pd,kd]=zpkdata(sys,'v')