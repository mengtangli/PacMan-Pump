% Plot moving 3D gear set 4 vertex with 2R-3R gears, no L
% Mengtang Li
% Apr 8
clear; clf; close all;

load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\x_record1.mat');
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\x_record2.mat');
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\y_record1.mat');
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\y_record2.mat');
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\z_record1.mat');
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\13\z_record2.mat');

figure('OuterPosition', [50 50 1000 1000]);
grid minor; grid on; hold on;
% plot x-axis
xaxis_x = 0:0.1:2;
xaxis_y = zeros(1,size(xaxis_x,2));
xaxis_z = zeros(1,size(xaxis_x,2));
% plot y-axis
yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.1:2;
yaxis_z = zeros(1,size(xaxis_x,2));
% plot z-axis
zaxis_x = zeros(1,size(xaxis_x,2));
zaxis_y = zeros(1,size(xaxis_x,2));
zaxis_z = 0:0.1:2;
% specify the lengh of axis
xlim([-3 3]); ylim([-3 3]); zlim([-3 3]);
plot3(xaxis_x,xaxis_y,xaxis_z,'k','linewidth',2);
grid on; grid minor; hold on;
plot3(yaxis_x,yaxis_y,yaxis_z,'k','linewidth',2);
plot3(zaxis_x,zaxis_y,zaxis_z,'k','linewidth',2);
% legend('x-axis', 'y-axis', 'z-axis');
text(2.0,0.05,0.05,'x');
text(0.05,2.0,0.05,'y');
text(0.05,0.05,2.0,'z');
daspect([1 1 1]);
view(24,24);

alpha = pi/6;
beta = pi/6;
R = 1;
% L = (2*sqrt(3)-3)*R;
h = 1;
angle = 0:0.02:4*pi;
n = size(angle,2);

plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle+sqrt(7)/2*R, 'r', 'linewidth',2)
plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle-sqrt(7)/2*R, 'r', 'linewidth',2)
plot3(x_record1, y_record1, z_record1, 'b-.','linewidth',1);
plot3(x_record2, y_record2, z_record2, 'r-.','linewidth',1);
pause(2);

for i = 1:1:n
    wt = angle(i);
    
    % cosa = 3/4; sina = sqrt(7)/4;
    Ainv = [ 3*cos(2*wt)/4,          3*sin(2*wt)/4,         7^(1/2)/4;
        -sin(2*wt),                cos(2*wt),             0;
        -7^(1/2)*cos(2*wt)/4,   -7^(1/2)*sin(2*wt)/4,       3/4;];
    
    Binv = [ cos(-1.5*wt)    sin(-1.5*wt)    0;
        -sin(-1.5*wt)    cos(-1.5*wt)    0;
        0           0        1;];
    
    C = Binv*Ainv;
    % Now the C is correct
    
    % distance between origin and center of rotating circle
    R_x = 0;
    R_y = 0;
    R_z = 0;
    
    % the rotating circle [blue]
    theta = 0:0.02:2*pi;
    rho_x = 2*R*(cos(theta)*C(1,1)+sin(theta)*C(2,1));
    rho_y = 2*R*(cos(theta)*C(1,2)+sin(theta)*C(2,2));
    rho_z = 2*R*(cos(theta)*C(1,3)+sin(theta)*C(2,3));
    % lower part half circle 1st point (0 deg) part
    phi = 0:0.02:pi/2;
    rho_p1_x = 2*R*cos(phi)*cos(0)*C(1,1)+2*R*cos(phi)*sin(0)*C(2,1)+2*R*sin(phi)*C(3,1);
    rho_p1_y = 2*R*cos(phi)*cos(0)*C(1,2)+2*R*cos(phi)*sin(0)*C(2,2)+2*R*sin(phi)*C(3,2);
    rho_p1_z = 2*R*cos(phi)*cos(0)*C(1,3)+2*R*cos(phi)*sin(0)*C(2,3)+2*R*sin(phi)*C(3,3);
    % lower part half circle 2nd point (90 deg) part
    rho_p2_x = 2*R*cos(phi)*cos(pi/2)*C(1,1)+2*R*cos(phi)*sin(pi/2)*C(2,1)+2*R*sin(phi)*C(3,1);
    rho_p2_y = 2*R*cos(phi)*cos(pi/2)*C(1,2)+2*R*cos(phi)*sin(pi/2)*C(2,2)+2*R*sin(phi)*C(3,2);
    rho_p2_z = 2*R*cos(phi)*cos(pi/2)*C(1,3)+2*R*cos(phi)*sin(pi/2)*C(2,3)+2*R*sin(phi)*C(3,3);
    % lower part half circle 3rd point (180 deg) pary
    rho_p3_x = 2*R*cos(phi)*cos(pi)*C(1,1)+2*R*cos(phi)*sin(pi)*C(2,1)+2*R*sin(phi)*C(3,1);
    rho_p3_y = 2*R*cos(phi)*cos(pi)*C(1,2)+2*R*cos(phi)*sin(pi)*C(2,2)+2*R*sin(phi)*C(3,2);
    rho_p3_z = 2*R*cos(phi)*cos(pi)*C(1,3)+2*R*cos(phi)*sin(pi)*C(2,3)+2*R*sin(phi)*C(3,3);
    % lower part half circle 4th point (270 deg) pary
    rho_p4_x = 2*R*cos(phi)*cos(3*pi/2)*C(1,1)+2*R*cos(phi)*sin(3*pi/2)*C(2,1)+2*R*sin(phi)*C(3,1);
    rho_p4_y = 2*R*cos(phi)*cos(3*pi/2)*C(1,2)+2*R*cos(phi)*sin(3*pi/2)*C(2,2)+2*R*sin(phi)*C(3,2);
    rho_p4_z = 2*R*cos(phi)*cos(3*pi/2)*C(1,3)+2*R*cos(phi)*sin(3*pi/2)*C(2,3)+2*R*sin(phi)*C(3,3);
    
    % vertex of normal vector of radius R in moving frame
    rho_xvt = 0+0+2*R*C(3,1);
    rho_yvt = 0+0+2*R*C(3,2);
    rho_zvt = 0+0+2*R*C(3,3);
    % 1st point (0 deg) in moving frame
    rho_x1 = 2*R*cos(0)*C(1,1)+2*R*sin(0)*C(2,1);
    rho_y1 = 2*R*cos(0)*C(1,2)+2*R*sin(0)*C(2,2);
    rho_z1 = 2*R*cos(0)*C(1,3)+2*R*sin(0)*C(2,3);
    % 2nd point (90 deg) in moving frame
    rho_x2 = 2*R*cos(pi/2)*C(1,1)+2*R*sin(pi/2)*C(2,1);
    rho_y2 = 2*R*cos(pi/2)*C(1,2)+2*R*sin(pi/2)*C(2,2);
    rho_z2 = 2*R*cos(pi/2)*C(1,3)+2*R*sin(pi/2)*C(2,3);
    % 3rd point (180 deg) in moving frame
    rho_x3 = 2*R*cos(pi)*C(1,1)+2*R*sin(pi)*C(2,1);
    rho_y3 = 2*R*cos(pi)*C(1,2)+2*R*sin(pi)*C(2,2);
    rho_z3 = 2*R*cos(pi)*C(1,3)+2*R*sin(pi)*C(2,3);
    % 4th point (270 deg) in moving frame
    rho_x4 = 2*R*cos(3*pi/2)*C(1,1)+2*R*sin(3*pi/2)*C(2,1);
    rho_y4 = 2*R*cos(3*pi/2)*C(1,2)+2*R*sin(3*pi/2)*C(2,2);
    rho_z4 = 2*R*cos(3*pi/2)*C(1,3)+2*R*sin(3*pi/2)*C(2,3);
    
    % center of moving circle in fixed frame
    r_x = R_x + rho_x;
    r_y = R_y + rho_y;
    r_z = R_z + rho_z;
    % vertex of normal vector of radius R in fixed frame
    r_xvt = R_x + rho_xvt;
    r_yvt = R_y + rho_yvt;
    r_zvt = R_z + rho_zvt;
    % 1st point (0 deg) in fixed frame
    r_x1 = R_x + rho_x1;
    r_y1 = R_y + rho_y1;
    r_z1 = R_z + rho_z1;
    % 2nd point (90 deg) in fixed frame
    r_x2 = R_x + rho_x2;
    r_y2 = R_y + rho_y2;
    r_z2 = R_z + rho_z2;
    % 3rd point (180 deg) in fixed frame
    r_x3 = R_x + rho_x3;
    r_y3 = R_y + rho_y3;
    r_z3 = R_z + rho_z3;
    % 4th point (270 deg) in fixed frame
    r_x4 = R_x + rho_x4;
    r_y4 = R_y + rho_y4;
    r_z4 = R_z + rho_z4;
    % lower part half circle 1st point (0 deg) part
    r_p1_x = R_x + rho_p1_x;
    r_p1_y = R_y + rho_p1_y;
    r_p1_z = R_z + rho_p1_z;
    % lower part half circle 2nd point (120 deg) part
    r_p2_x = R_x + rho_p2_x;
    r_p2_y = R_y + rho_p2_y;
    r_p2_z = R_z + rho_p2_z;
    % lower part half circle 3rd point (240 deg) part
    r_p3_x = R_x + rho_p3_x;
    r_p3_y = R_y + rho_p3_y;
    r_p3_z = R_z + rho_p3_z;
    % lower part half circle 4th point (270 deg) part
    r_p4_x = R_x + rho_p4_x;
    r_p4_y = R_y + rho_p4_y;
    r_p4_z = R_z + rho_p4_z;
    
    plotcircle1 = plot3(r_x,r_y,r_z,'b','linewidth',2);
    plotcircle2 = plot3(r_p1_x,r_p1_y,r_p1_z,'r','linewidth',2);
    plotcircle3 = plot3(r_p2_x,r_p2_y,r_p2_z,'r','linewidth',2);
    plotcircle4 = plot3(r_p3_x,r_p3_y,r_p3_z,'r','linewidth',2);
    plotcircle5 = plot3(r_p4_x,r_p4_y,r_p4_z,'r','linewidth',2);
    plotline0 = plot3([R_x; r_xvt], [R_y; r_yvt], [R_z; r_zvt],'b','linewidth',2);
    plotline1 = plot3([r_x1; r_x3], [r_y1; r_y3], [r_z1; r_z3],'r','linewidth',2);
    plotline2 = plot3([r_x2; r_x4], [r_y2; r_y4], [r_z2; r_z4],'r','linewidth',2);
            
    pause(0.05);
    if (i < n)
        delete(plotcircle1);
        delete(plotcircle2);
        delete(plotcircle3);
        delete(plotcircle4);
        delete(plotcircle5);
        delete(plotline0);
        delete(plotline1);
        delete(plotline2);   
    end
    
end