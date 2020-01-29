function [worldOrientation, worldLocation, H, A] = PoseEstimation(worldPoints, imagePoints, K)
%Let's fit the worldPoints with the imagePoints
H=imagePoints/worldPoints;

%partial_sum=0;
% for i = 1:size(imagePoints,2)
%     partial_sum=sum([abs(imagePoints(:,i) - H*worldPoints(:,i))]./imagePoints(:,i)) + partial_sum
% end
% 
% error=partial_sum/(2*size(imagePoints,2))

KO=[K zeros(3,1)];
A=KO\H;


worldOrientation=A(1:3,1:3);
worldLocation=A(1:3,4);

end

