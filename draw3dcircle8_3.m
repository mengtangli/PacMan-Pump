% Plot moving 3D gear set with 1.5R-2R gears, no L
% Mengtang Li
% Mar 30

% OK this is wrong. Check version 9

clear; clf; close all;

figure('OuterPosition', [50 50 1000 1000]);
grid minor; grid on; hold on;
% plot axes
xaxis_x = 0:0.1:2;
xaxis_y = zeros(1,size(xaxis_x,2));
xaxis_z = zeros(1,size(xaxis_x,2));

yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.1:2;
yaxis_z = zeros(1,size(xaxis_x,2));

zaxis_x = zeros(1,size(xaxis_x,2));
zaxis_y = zeros(1,size(xaxis_x,2));
zaxis_z = 0:0.1:2;

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
angle = 0:0.01:2*pi;
n = size(angle,2);
x_record1 = zeros(1,n);
y_record1 = zeros(1,n);
z_record1 = zeros(1,n);
x_record2 = zeros(1,n);
y_record2 = zeros(1,n);
z_record2 = zeros(1,n);
x_record3 = zeros(1,n);
y_record3 = zeros(1,n);
z_record3 = zeros(1,n);

% plot3(L/2*cos(angle), L/2*sin(angle), 0*angle, 'r', 'linewidth',2)
plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle+sqrt(7)/2*R, 'r', 'linewidth',2)
plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle-sqrt(7)/2*R, 'r', 'linewidth',2)
pause(2);

for i = 1:1:n
    wt = angle(i);
    % cosa = 3/4; sina = sqrt(7)/4; 
    Ainv = [ (3*cos(2*wt))/4,        (3*sin(2*wt))/4,         7^(1/2)/4;
             -sin(2*wt),                cos(2*wt),               0;
         -(7^(1/2)*cos(2*wt))/4,    -(7^(1/2)*sin(2*wt))/4,     3/4];
    
    Binv = [ cos(-1.5*wt) sin(-1.5*wt) 0;
        -sin(-1.5*wt) cos(-1.5*wt) 0;
        0       0       1;];
    
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
    
    % vertex of normal vector of radius R in moving frame
    rho_xvt = 0+0+2*R*C(3,1);
    rho_yvt = 0+0+2*R*C(3,2);
    rho_zvt = 0+0+2*R*C(3,3);
    % 1st point (0 deg) in moving frame
    rho_x1 = 2*R*(cos(0)*C(1,1)+sin(0)*C(2,1));
    rho_y1 = 2*R*(cos(0)*C(1,2)+sin(0)*C(2,2));
    rho_z1 = 2*R*(cos(0)*C(1,3)+sin(0)*C(2,3));
    % 2nd point (90 deg) in moving frame
    rho_x2 = 2*R*(cos(pi/2)*C(1,1)+sin(pi/2)*C(2,1));
    rho_y2 = 2*R*(cos(pi/2)*C(1,2)+sin(pi/2)*C(2,2));
    rho_z2 = 2*R*(cos(pi/2)*C(1,3)+sin(pi/2)*C(2,3));
    % 3rd point (180 deg) in moving frame
    rho_x3 = 2*R*(cos(pi)*C(1,1)+sin(pi)*C(2,1));
    rho_y3 = 2*R*(cos(pi)*C(1,2)+sin(pi)*C(2,2));
    rho_z3 = 2*R*(cos(pi)*C(1,3)+sin(pi)*C(2,3));
    % 4th point (270 deg) in moving frame
    rho_x4 = 2*R*(cos(3*pi/2)*C(1,1)+sin(3*pi/2)*C(2,1));
    rho_y4 = 2*R*(cos(3*pi/2)*C(1,2)+sin(3*pi/2)*C(2,2));
    rho_z4 = 2*R*(cos(3*pi/2)*C(1,3)+sin(3*pi/2)*C(2,3));
    
    % center of moving circle in fixed frame
    r_x = R_x + rho_x;
    r_y = R_y + rho_y;
    r_z = R_z + rho_z;
    %vertex of normal vector of radius R in fixed frame
    r_xvt = R_x + rho_xvt;  x_record1(i) = r_xvt;
    r_yvt = R_y + rho_yvt;  y_record1(i) = r_yvt;
    r_zvt = R_z + rho_zvt;  z_record1(i) = r_zvt;
    % 1st point (0 deg) in fixed frame
    r_x1 = R_x + rho_x1;  x_record2(i) = r_x1;
    r_y1 = R_y + rho_y1;  y_record2(i) = r_y1;
    r_z1 = R_z + rho_z1;  z_record2(i) = r_z1;
    % 2nd point (90 deg) in fixed frame
    r_x2 = R_x + rho_x2;  x_record3(i) = r_x2;
    r_y2 = R_y + rho_y2;  y_record3(i) = r_y2;
    r_z2 = R_z + rho_z2;  z_record3(i) = r_z2;
    % 3rd point (180 deg) in fixed frame
    r_x3 = R_x + rho_x3;
    r_y3 = R_y + rho_y3;
    r_z3 = R_z + rho_z3;
    % 4th point (270 deg) in fixed frame
    r_x4 = R_x + rho_x4;
    r_y4 = R_y + rho_y4;
    r_z4 = R_z + rho_z4;
    
    plotcircle1 = plot3(r_x,r_y,r_z,'b','linewidth',2);
    plotline0 = plot3([R_x; r_xvt], [R_y; r_yvt], [R_z; r_zvt],'b','linewidth',2);
    plotline1 = plot3([r_x1; r_x3], [r_y1; r_y3], [r_z1; r_z3],'y','linewidth',2);
    plotline2 = plot3([r_x2; r_x4], [r_y2; r_y4], [r_z2; r_z4],'g','linewidth',2);
    
    
    pause(0.02);
    if (i < n)
        delete(plotcircle1);
        delete(plotline0);
        delete(plotline1);
        delete(plotline2);
    end
    
end
