function solution = IntersectionLineConic(C1,l)
a1= C1(1,1);
b1=2*C1(1,2);
c1=C1(2,2);
d1=2*C1(1,3);
e1=2*C1(2,3);
f1=C1(3,3);

syms x y;
eq1 = a1*x^2 + b1*x*y + c1*y^2 + d1*x + e1*y + f1; 
eq2 = l(1)*x+l(2)*y+l(3);
eqs=[eq1 eq2];

solution = solve(eqs, [x y]);
solution = [double(solution.x'); double(solution.y'); ones(1,2)];


end

