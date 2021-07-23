function E = Efunc(I1,I2,I3,J,K,Mb,Mw,d,dq1,dq2,dq3,dq4,g,l,q3,q4,r)
%EFUNC
%    E = EFUNC(I1,I2,I3,J,K,MB,MW,D,DQ1,DQ2,DQ3,DQ4,G,L,Q3,Q4,R)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    05-Sep-2019 14:17:35

t2 = sin(q3);
t7 = cos(q4);
t8 = sin(q4);
t18 = t2.^2;
t3 = dq1.*t7+dq2.*t8+dq4.*l.*t18;
t4 = l.^2;
t5 = dq3.^2;
t6 = dq4.^2;
t10 = d.*dq4;
t11 = dq1.*t7.*2.0;
t12 = dq2.*t8.*2.0;
t9 = t10+t11+t12;
t13 = -t10+t11+t12;
t14 = cos(q3);
t15 = t9.^2;
t16 = 1.0./r.^2;
t17 = t13.^2;
E = I2.*t5.*(1.0./2.0)+K.*t6+Mw.*t15.*(1.0./8.0)+Mw.*t17.*(1.0./8.0)+Mb.*(t4.*t5+t3.^2-t4.*t6.*(cos(q3.*4.0).*(1.0./8.0)-1.0./8.0)).*(1.0./2.0)+I1.*t6.*t18.*(1.0./2.0)+J.*t15.*t16.*(1.0./8.0)+J.*t16.*t17.*(1.0./8.0)+I3.*t6.*t14.^2.*(1.0./2.0)+Mb.*g.*l.*t14;
