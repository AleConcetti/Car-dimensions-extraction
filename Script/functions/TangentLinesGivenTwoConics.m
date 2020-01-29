function solution = TangentLinesGivenTwoConics(C1,C2)
C1star=inv(C1);
C2star=inv(C2);

a1= C1star(1,1);
b1=2*C1star(1,2);
c1=C1star(2,2);
d1=2*C1star(1,3);
e1=2*C1star(2,3);
f1=C1star(3,3);

a2= C2star(1,1);
b2=2*C2star(1,2);
c2=C2star(2,2);
d2=2*C2star(1,3);
e2=2*C2star(2,3);
f2=C2star(3,3);

syms x y;
eq1= a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1; 
eq2= a2*x^2 + b2*x*y + c2*y^2 + d2*x + e2*y + f2;
eqs=[eq1 eq2];

solution = solve(eqs, [x y]);
solution = [double(solution.x'); double(solution.y'); ones(1,4)];

end

