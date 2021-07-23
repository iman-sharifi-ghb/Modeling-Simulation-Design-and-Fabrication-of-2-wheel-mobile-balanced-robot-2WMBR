clc
clear 
close all
[A,B,C,D] = linearize_2wbr(12);
sys  = ss(A,B,C,D);
sys = tf(sys)