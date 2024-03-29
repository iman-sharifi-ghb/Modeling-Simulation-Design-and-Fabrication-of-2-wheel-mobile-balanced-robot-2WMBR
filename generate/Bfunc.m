function B = Bfunc(I1,I3,J,Mb,Mw,dq1,dq2,dq4,g,l,q3,q4,r)
%BFUNC
%    B = BFUNC(I1,I3,J,MB,MW,DQ1,DQ2,DQ4,G,L,Q3,Q4,R)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    05-Sep-2019 14:17:35

t2 = dq4.^2;
t3 = cos(q3);
t4 = cos(q4);
t5 = sin(q4);
t6 = r.^2;
t7 = sin(q3);
B = [0.0;0.0;-t7.*(Mb.*g.*l+I1.*t2.*t3-I3.*t2.*t3+Mb.*l.^2.*t2.*t3+Mb.*dq1.*dq4.*l.*t3.*t4.*2.0+Mb.*dq2.*dq4.*l.*t3.*t5.*2.0);1.0./r.^2.*(dq1.*t5-dq2.*t4).*(J.*dq1.*t4.*2.0+J.*dq2.*t5.*2.0+Mb.*dq1.*t4.*t6+Mb.*dq2.*t5.*t6+Mw.*dq1.*t4.*t6.*2.0+Mw.*dq2.*t5.*t6.*2.0+Mb.*dq4.*l.*t6.*t7.^2);0.0;0.0];
