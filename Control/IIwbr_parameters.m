function p = IIwbr_parameters()
% 2wbr parameters---------------------------------------------------
p.d = 0.22;
p.l = 0.048;
p.r = 0.05;
p.Mb = 1.021;
p.Mw = 0.08;
p.g = 9.81;
p.I1 = 5.493991e-3;
p.I2 = 2.106395e-3;
p.I3 = 4.001042e-3;
p.J = 1.14927e-4; 
p.K = 6.9505e-5;
p.Ca = 0.05;
% dc motor parameters--------------------------------------------------
p.Km = 1.34;
p.Kb = 1;
p.R = 1.7;
p.L = 0.05;
J_motor = 4.5e-6;  
J_gearbox = 1.5842e-5;
N = 59;
p.J_e = J_gearbox+N^2*J_motor;
p.B = 0.043;
end