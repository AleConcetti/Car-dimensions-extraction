function C1 = FindConicFrom5Points(x1,y1)
A1=[x1.^2 x1.*y1 y1.^2 x1 y1 ones(size(x1))];
N1 = null(A1);
[a b c d e f]=deal(N1(1),N1(2),N1(3),N1(4),N1(5),N1(6));
C1=[a b/2 d/2; b/2 c e/2; d/2 e/2 f];
end

