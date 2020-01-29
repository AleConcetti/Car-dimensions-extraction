function mid = MiddlePointByCR(p1,p2,v)
%This function return the image of the real middle point in the image
%plane, by putting the cross ratio equall to -1.

% v is the vanishing point
% m is the middle point
% p1, p2 are the symmetric points
line=cross(p1,p2);
line=line/line(3);

syms m1 m2;
m=[m1;m2;1];

m_belong_to_line = line'*m;
CR  = (-Distance(m,p1)/Distance(m,p2))/(Distance(v,p1)/Distance(v,p2)) + 1;
[mid1 mid2] = solve([CR, m_belong_to_line], [m1 m2]);

if mid1(1)>0 && mid2(1)>0
    mid1=mid1(1);
    mid2=mid2(1);
elseif mid1(2)>0 && mid2(2)>0
    mid1=mid1(2);
    mid2=mid2(2);
else
    error('Error: the middle point found is not in the image') 
    mid1= -1;
    mid2= -1;
end

mid=[double(mid1) ; double(mid2) ; 1];

end

