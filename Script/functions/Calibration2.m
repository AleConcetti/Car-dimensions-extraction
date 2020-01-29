function K = Calibration2(I,J,vh,vv,vl)
syms a b c d e f;
w=[ a b d;...
    b c e;...
    d e f];

eq1=I'*w*I;
eq2=J'*w*J;
eq3=vh'*w*vv;
eq4=vv'*w*vl;
eq5=vl'*w*vh;
eq6=f-1;

eqs=[eq1;eq2;eq3;eq4;eq5;eq6];

par=solve(eqs, [a b c d e f]);

w=[ double(par.a) double(par.b) double(par.d);...
    double(par.b) double(par.c) double(par.e);...
    double(par.d) double(par.e) double(par.f) ];

w
inv(w)
eig(inv(w))
K=chol(inv(w),'upper');

end
