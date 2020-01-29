%% SETUP
close all
clear all
clc

addpath(genpath('functions'));
addpath(genpath('images'));
addpath(genpath('data'));

fprintf("I'm setting up the image\n")
im=SetupImage(imread('Image1.jpg'));

% RESCALING THE PICTURE FOR FASTER CALCULATION
sizeIm = size(im);
lenIm= sizeIm(1);
fprintf("Set-up completed\n\n")
%% FEATURES SELECTION (ellipses and symmetric points)
fprintf("I'm searching for the ellipes.\n")
[C_ant, C_post] = FindWheelEllipses(im);
C1=C_ant;
C2=C_post;
fprintf("Ellipses found! :)\n\n")

fprintf("I'm searching for symmetric points.\n")
symPoints=SymPointsSelection(im);

fprintf("I've found some pairs of symmetric points! :) \n\n")

figure(1), imshow(im,[]), hold on,
plot(symPoints(:,1), symPoints(:,2),'g+', 'LineWidth', 4)
title('Symmetric features')

%PLOT THE FIGURE
conic1=BuildMatrixFrom('conic',C1,lenIm);
conic2=BuildMatrixFrom('conic',C2,lenIm);
figure(2),
imshow(conic1+conic2+im, [0 1]);

%% WHEELS' PLANE RECTIFICATION

%FIND THE LINES TANGENT TO THE TWO CONICS FOUND
solution = TangentLinesGivenTwoConics(C1,C2);
l1=solution(:,1);
l2=solution(:,2);
l3=solution(:,3);
l4=solution(:,4);

%PLOT THEM
line1=BuildMatrixFrom('line',l1,lenIm);
line2=BuildMatrixFrom('line',l2,lenIm);
line3=BuildMatrixFrom('line',l3,lenIm);
line4=BuildMatrixFrom('line',l4,lenIm);
sol2=conic1+conic2+im+line1+line4;
sol1=conic1+conic2+im+line2+line3;
figure(3),
imshow([sol1 sol2], [0 1]), title('I optain 2 pair of lines');

% FIND TANGENT POINTS
left1=IntersectionLineConic(C1,l2);
left1=left1(:,1); %ruota posteriore up
left2=IntersectionLineConic(C1,l3);
left2=left2(:,1); %ruota posteriore down
right1=IntersectionLineConic(C2,l2);
right1=right1(:,1); %ruota anteriore up
right2=IntersectionLineConic(C2,l3);
right2=right2(:,1); %ruota anteriore down
intersectionPoints = BuildMatrixFrom('point', left1,lenIm)+ BuildMatrixFrom('point', left2,lenIm)+ BuildMatrixFrom('point', right1,lenIm)+ BuildMatrixFrom('point', right2,lenIm);
%this matrix below is used oly for plotting
hor_vv_plot = conic1+conic2+im+line2+line3+intersectionPoints;

% FIND CIRCULAR POINTS
%Vertical Vanishing point
line_left=cross(left1,left2);
line_left=line_left/line_left(3);
line_right=cross(right1,right2);
line_right=line_right/line_right(3);

ll=BuildMatrixFrom('line',line_left,lenIm);
lr=BuildMatrixFrom('line',line_right,lenIm);

vet_vv_plot= conic1+conic2+im+intersectionPoints+ll+lr;
figure(4),
imshow([vet_vv_plot,hor_vv_plot], [0 1]);

vv=cross(line_left,line_right);
vv=vv/vv(3);
%Horizontal Vanishing point
vh=cross(l2,l3);
vh=vh/vh(3);


%Find the line at infinity as the intersection of vanishing points
l_inf=cross(vv,vh);
l_inf=l_inf/l_inf(3);
%Find the circular points as the intersection between an ellipse and the
%l_inf
circPoints = IntersectionLineConic(C2,l_inf);
I=circPoints(:,1);
J=circPoints(:,2);

%Find Hr with the svd function
dualCinf=I*J'+J*I';
[U S V]=svd(dualCinf)
S1=[(S(1,1))^(0.5)         0             0; ...
           0         (S(2,2))^(0.5)      0; ...
           0               0             1    ];
% Let dCinf=[1 0 0; 0 1 0; 0 0 0]
% Note that S = S1*dCinf*S1, hence...
% dualCinf=(U*S1)*dCinf*(V*S1)'
Hr=inv(U*S1);
Hr=Hr/Hr(3,3);

% TROVA DIAMETRO RUOTA POSTERIORE
diam_left = [left1 left2];

%TROVA DIAMETRO RUOTA ANTERIORE
diam_right = [right1 right2];

%INTERSECA I DUE PUNTI MEDI DEI DIAMETRI
med_left = MiddlePointByCR(diam_left(:,1),diam_left(:,2),vh);
med_left=med_left/med_left(3);
med_right = MiddlePointByCR(diam_right(:,1),diam_right(:,2),vh);
med_right=med_right/med_right(3);

%Now let's transform all the points and found the ,asked ratio
%Transform the wheel-to-wheel distance
med_left_trasf = Hr*med_left;
med_left_trasf=med_left_trasf/med_left_trasf(3);

med_right_trasf = Hr*med_right;
med_right_trasf=med_right_trasf/med_right_trasf(3);

lenCar = [med_left_trasf med_right_trasf];

%Transform the diameter of wheel circle
left1_trasf=Hr*left1;
left1_trasf=left1_trasf/left1_trasf(3);

left2_trasf=Hr*left2;
left2_trasf=left2_trasf/left2_trasf(3);

diamLtrasf = [left1_trasf left2_trasf];


right1_trasf=Hr*right1;
right1_trasf=right1_trasf/right1_trasf(3);

right2_trasf=Hr*right2;
right2_trasf=right2_trasf/right2_trasf(3);

diamRtrasf = [right1_trasf right2_trasf];

fprintf('CHECK: The diameter of the two wheel circle are the same after transformation')
Lenght(diamLtrasf)
Lenght(diamRtrasf)

%Make the ratio!
fprintf('---->THE RATIO BETWEEN THE DIAMETER OF WHEEL AND WHEEL-TO-WHEEL DISTANCE IS:')
result= Lenght(diamLtrasf)/Lenght(lenCar)


%% ------------------------- CALIBRATE THE CAMERA--------------------------
%The camera is zero-screw but not natural so I need at least
%4 equations for finding the K matrix(fx,fy,Ux,Uy)

%I use the vanisching points to find 4 equations,

%First of all I need the lateral vanishing points. Let's use some symmetric features:
vl = LateralVanishingPoints(symPoints);

%Now I need another pair of orthogonal vaishing points. Let's consider the
%rear rim:
%If a circonference is crossed by a vertical line and an horizontal one, both
%passing through the center, then the intersection points are 4 (on the top, on the bottom, on the right and on the left)
% Interpoling the point as follows I can find a pair of ortogonal lines
line_joining_centers = cross(med_left, med_right);
line_joining_centers = Normalize("vector", line_joining_centers);

sol = IntersectionLineConic(C2, line_joining_centers);

p_left = sol(:,2);
p_right = sol(:,1);
p_down = right1
p_up = right2


v1=cross(cross(p_up,p_right),cross(p_down,p_left));
v1=Normalize("vector", v1);

v2=cross(cross(p_up,p_left),cross(p_down,p_right));
v2=Normalize("vector", v2);

%Now I have 5 vanishing points:
%three of them are orthogonal each other: vh, vv, vl (horizontal, vertical, lateral)
%So the 3 equations are:
%eq1 = vh'*w*vv;
%eq2 = vh'*w*vl;
%eq3 = vl'*w*vv;
%Then I have v1 and v2 that are orthogonal, so:
%eq4 = v1'*w*v2

[K iac] = Calibration(vh,vv,vl,v1,v2)

%% ---------------------RECONSTRUCT 3D POINTS RECONSTRUCTION--------------------
% Fix the reference frame at a suitable position on the
% symmetry plane of the car: I choose to put it in the middle point of the
% lower pair of symmetric points p1 , p2

p1=[symPoints(5,:)'; 1]; % on the left
p2=[symPoints(6,:)'; 1]; % on the right
o=MiddlePointByCR(p1,p2, vl);
% o is my reference frame center
% (x // short side; y // long side; z is the height)
% terna destrorsa con z verso l'alto
op=BuildMatrixFrom("point",o,lenIm);
figure(6), imshow(im+op), title("Center of the reference frame");

%  - Consider 2 planes:
%       - plane_h parallel to the floor
%       - plane-v with axis parallel to the long side of the car
%  - For each plane let's find the circular points (since we know iac from calibration):

%PLANE H
line_inf=cross(vl,vh);
line_inf=line_inf/line_inf(3);
circPointsH = IntersectionLineConic(iac,line_inf);
I=circPointsH(:,1);
J=circPointsH(:,2);
dualCinf=I*J'+J*I';
[U S V]=svd(dualCinf);
S1=[(S(1,1))^(0.5)         0             0; ...
           0         (S(2,2))^(0.5)      0; ...
           0               0             1    ];
Hh=inv(U*S1);
Hh=Hh/Hh(3,3);

%PLANE V
line_inf=cross(vl,vv);
line_inf=line_inf/line_inf(3);
circPointsV = IntersectionLineConic(iac,line_inf);
I=circPointsV(:,1);
J=circPointsV(:,2);
dualCinf=I*J'+J*I';
[U S V]=svd(dualCinf);
S1=[(S(1,1))^(0.5)         0             0; ...
           0         (S(2,2))^(0.5)      0; ...
           0               0             1    ];
Hv=inv(U*S1);
Hv=Hv/Hv(3,3);

%Now I scale the transformation matrix Hv in such a way Hv and Hh transform
%the segment p1-p2 in two segment with the same lenght.

segv=Lenght(Normalize("segment", Hv*[p1, p2]))
segh=Lenght(Normalize("segment", Hh*[p1, p2]))

Hv=[segh/segv 0 0; 0 segh/segv 0; 0 0 1]*Hv;

segv=Lenght(Normalize("segment", Hv*[p1, p2]))
segh=Lenght(Normalize("segment", Hh*[p1, p2]))

worldPoints=[];
n=size(symPoints,1);
for i= 1:2:n
    % i= 1(I), 3(II), 5(III), 7(IV), 9(V)
    % Symmetric points selection
    p1 = [symPoints(i,:)'; 1]; % on the left
    p2 = [symPoints(i+1,:)'; 1]; % on the right
    m = MiddlePointByCR(p1,p2,vl);
    x1 = [m,p1];
    x2 = [m,p2];

    vert_axis_through_m=cross(m,vv);
    vert_axis_through_m=vert_axis_through_m/vert_axis_through_m(3);

    y_axis=cross(o,vh);
    y_axis=y_axis/y_axis(3);

    h = cross( y_axis, vert_axis_through_m );
    h=h/h(3);

    z1 = [m,h];
    z2 = z1;

    y1 = [h,o];
    y2 = y1;


    p=cross(cross(o,vv),cross(m,vh));
    p=Normalize("vector", p);
    s=Distance(vh,p)/Distance(vh,m); %SCALE FACTOR DUE TO PROSPECTIVE
    scaling=[s 0 0; 0 s 0; 0 0 1];
    %Let's transform the segments
    segment_x1=Normalize("segment", scaling*Hv*x1);
    segment_y1=Normalize("segment", Hh*y1);
    segment_z1=Normalize("segment", scaling*Hv*z1);

    segment_x2=Normalize("segment", scaling*Hv*x2);
    segment_y2=Normalize("segment", Hh*y2);
    segment_z2=Normalize("segment", scaling*Hv*z2);

    %Let's find the real coordinates

    realx1=  -Distance(segment_x1(:,1),segment_x1(:,2));
    realy1=   Distance(segment_y1(:,1),segment_y1(:,2));
    realz1=   Distance(segment_z1(:,1),segment_z1(:,2));

    p_left=[realx1; realy1; realz1; 1];

    realx2=   Distance(segment_x2(:,1),segment_x2(:,2));
    realy2=   Distance(segment_y2(:,1),segment_y2(:,2));
    realz2=   Distance(segment_z2(:,1),segment_z2(:,2));

    p_right=[realx2; realy2; realz2; 1];

    worldPoints = [worldPoints, p_left, p_right];

end

fprintf('---->THE COORDINATES OF THE SYMMETRIC POINTS ARE THE FOLLOWINGS:')
worldPoints

pcshow((worldPoints(1:3,:))','VerticalAxis','Y','VerticalAxisDir','down', ...
 'MarkerSize',1000);
xlabel('X');
ylabel('Y');
zlabel('Z');

imagePoints=[symPoints'; ones(1, size(symPoints,1))];

%% -------------------- CAMERA POSE ESTIMATION ---------------------------


[worldOrientation, worldLocation, H, A] = PoseEstimation(worldPoints,imagePoints,K);

Ainv=[worldOrientation, -worldOrientation'*worldLocation];
Ainv=[Ainv; 0 0 0 1]

camera=Ainv*[0;0;0;1]

