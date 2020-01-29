function vl = LateralVanishingPoints(symPoints)
%This function evaluate the lateral vanishing point making a mean of all
%the vanishing points found taking two generic pairs of symmetric points.
    n=size(symPoints,1);
    setOfVl=[];
    for i=1:2:n
        for j=i+2:2:n
            symPoints1=[symPoints(i,:)',symPoints(i+1,:)';ones(1,2)];
            symPoints2=[symPoints(j,:)',symPoints(j+1,:)';ones(1,2)];
            line1 = cross(symPoints1(:,1),symPoints1(:,2));
            line1=line1/line1(3);
            line2 = cross(symPoints2(:,1),symPoints2(:,2));
            line2=line2/line2(3);
            v=cross(line1,line2);
            v=v/v(3);
            setOfVl=[setOfVl , v];    
        end
    end
    
    m=size(setOfVl,2);
    vl = [  sum(setOfVl(1,:))/m;...
            sum(setOfVl(2,:))/m;...
                1                ];

end

