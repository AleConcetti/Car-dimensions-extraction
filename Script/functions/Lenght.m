function len = Lenght(segment)
p1=segment(:,1);
p2=segment(:,2);
len = ((p1(1)-p2(1))^2+(p1(2)-p2(2))^2)^0.5;
end

