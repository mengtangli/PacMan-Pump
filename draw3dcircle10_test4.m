% Plot moving 3D gear set with same size but lower part 
% (show lower semi circle). The center of the smaller half piece
% is rotating.
% Mengtang Li
% Apr 13
clear; clf; close all;

load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\10\x_re1.mat')
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\10\y_re1.mat')
load('C:\Users\lim14\Documents\MATLAB\Project_504\TAH\draw3dcircle\10\z_re1.mat')

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

plot3(x_re1, y_re1, z_re1, 'm-.','linewidth',1);

alpha = pi/6;
beta = pi/6;
R = 1;
L = (sqrt(3)-1)*2*R;
h = L;
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
x_record4 = zeros(1,n);
y_record4 = zeros(1,n);
z_record4 = zeros(1,n);
x_record5 = zeros(1,n);
y_record5 = zeros(1,n);
z_record5 = zeros(1,n);
x_record6 = zeros(1,n);
y_record6 = zeros(1,n);
z_record6 = zeros(1,n);

plot3(L/2*cos(angle), L/2*sin(angle), 0*angle, 'r', 'linewidth',2)
plot3(R*cos(angle), R*sin(angle), 0*angle+R, 'r', 'linewidth',2)
plot3(0.3*R*sin(pi/6)*cos(angle), 0.3*R*sin(pi/6)*sin(angle), -0.3*R*cos(pi/6)+0*angle, 'm-.','linewidth',1);

pause(2);

for i = 1:1:n
    wt = angle(i);
    
    Ainv = [ (3^(1/2)*cos(2*wt))/2,  (3^(1/2)*sin(2*wt))/2,    1/2;
        -sin(2*wt),               cos(2*wt),              0;
        -cos(2*wt)/2,            -sin(2*wt)/2,         3^(1/2)/2;];
    
    Binv = [ cos(-wt) sin(-wt) 0;
        -sin(-wt) cos(-wt) 0;
        0       0       1;];
    
    C = Binv*Ainv;
    % Now the C is correct
    
    % distance between origin and center of rotating circle
    R_x = -L/2*cos(2*wt);
    R_y = -L/2*sin(2*wt);
    R_z = 0;
    
    % the rotating circle [blue]
    theta = 0:0.02:2*pi;
    rho_x = 2*R*(cos(theta)*C(1,1)+sin(theta)*C(2,1));
    rho_y = 2*R*(cos(theta)*C(1,2)+sin(theta)*C(2,2));
    rho_z = 2*R*(cos(theta)*C(1,3)+sin(theta)*C(2,3));
    
    % lower part half circle [green]
    phi = 0:0.02:pi;
    rho2_x = 2*R*(cos(phi)*C(2,1)+sin(phi)*C(3,1));
    rho2_y = 2*R*(cos(phi)*C(2,2)+sin(phi)*C(3,2));
    rho2_z = 2*R*(cos(phi)*C(2,3)+sin(phi)*C(3,3));
    
    % bigger lower part half circle [cayon]
    phi = pi:0.02:2*pi;
    rho2_xe = 3*R*cos(phi)*C(2,1)+(3*R*sin(phi)-h)*C(3,1);
    rho2_ye = 3*R*cos(phi)*C(2,2)+(3*R*sin(phi)-h)*C(3,2);
    rho2_ze = 3*R*cos(phi)*C(2,3)+(3*R*sin(phi)-h)*C(3,3);
    
    % lowered bigger lower part half circle [cayon]
    phi2 = (pi+asin(0.1)):0.02:(2*pi-asin(0.1));
    rho2_xe2 = 3*R*cos(phi2)*C(2,1)+(3*R*sin(phi2)-h)*C(3,1);
    rho2_ye2 = 3*R*cos(phi2)*C(2,2)+(3*R*sin(phi2)-h)*C(3,2);
    rho2_ze2 = 3*R*cos(phi2)*C(2,3)+(3*R*sin(phi2)-h)*C(3,3);
    
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
    % extended 2nd point (90 deg) in moving frame
    rho_x2e = 3*R*cos(pi/2)*C(1,1)+3*R*sin(pi/2)*C(2,1)-h*C(3,1);
    rho_y2e = 3*R*cos(pi/2)*C(1,2)+3*R*sin(pi/2)*C(2,2)-h*C(3,2);
    rho_z2e = 3*R*cos(pi/2)*C(1,3)+3*R*sin(pi/2)*C(2,3)-h*C(3,3);
    % lowered extended 2nd point (90 deg) in moving frame
    rho_x2e2 = 3*R*cos(pi+asin(0.1))*C(2,1)+(3*R*sin(pi+asin(0.1))-h)*C(3,1);
    rho_y2e2 = 3*R*cos(pi+asin(0.1))*C(2,2)+(3*R*sin(pi+asin(0.1))-h)*C(3,2);
    rho_z2e2 = 3*R*cos(pi+asin(0.1))*C(2,3)+(3*R*sin(pi+asin(0.1))-h)*C(3,3);
    % extended 4nd point (270 deg) in moving frame
    rho_x4e = 3*R*cos(3*pi/2)*C(1,1)+3*R*sin(3*pi/2)*C(2,1)-h*C(3,1);
    rho_y4e = 3*R*cos(3*pi/2)*C(1,2)+3*R*sin(3*pi/2)*C(2,2)-h*C(3,2);
    rho_z4e = 3*R*cos(3*pi/2)*C(1,3)+3*R*sin(3*pi/2)*C(2,3)-h*C(3,3);
    % lowered extended 4nd point (270 deg) in moving frame
    rho_x4e2 = 3*R*cos(2*pi-asin(0.1))*C(2,1)+(3*R*sin(2*pi-asin(0.1))-h)*C(3,1);
    rho_y4e2 = 3*R*cos(2*pi-asin(0.1))*C(2,2)+(3*R*sin(2*pi-asin(0.1))-h)*C(3,2);
    rho_z4e2 = 3*R*cos(2*pi-asin(0.1))*C(2,3)+(3*R*sin(2*pi-asin(0.1))-h)*C(3,3);
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
    % lower part half circle
    r2_x = R_x + rho2_x;
    r2_y = R_y + rho2_y;
    r2_z = R_z + rho2_z;
    % bigger lower part half circle
    r2_xe = R_x + rho2_xe;
    r2_ye = R_y + rho2_ye;
    r2_ze = R_z + rho2_ze + (3^0.5)/2*h;
    % lowered bigger lower part half circle
    r2_xe2 = R_x + rho2_xe2;
    r2_ye2 = R_y + rho2_ye2;
    r2_ze2 = R_z + rho2_ze2 + (3^0.5)/2*h;    
    % extended 2nd point (90 deg) in moving frame
    r_x2e = R_x + rho_x2e;  x_record4(i) = r_x2e;
    r_y2e = R_y + rho_y2e;  y_record4(i) = r_y2e;
    r_z2e = R_z + rho_z2e + (3^0.5)/2*h;  z_record4(i) = r_z2e;
    % lowered extended 2nd point (90 deg) in moving frame
    r_x2e2 = R_x + rho_x2e2;  x_record5(i) = r_x2e2;
    r_y2e2 = R_y + rho_y2e2;  y_record5(i) = r_y2e2;
    r_z2e2 = R_z + rho_z2e2 + (3^0.5)/2*h;  z_record5(i) = r_z2e2;    
    % extended 4th point (270 deg) in moving frame
    r_x4e = R_x + rho_x4e;  %x_record5(i) = r_x4e;
    r_y4e = R_y + rho_y4e;  %y_record5(i) = r_y4e;
    r_z4e = R_z + rho_z4e + (3^0.5)/2*h;  %z_record5(i) = r_z4e;
    % lowered extended 4th point (270 deg) in moving frame
    r_x4e2 = R_x + rho_x4e2;  %x_record5(i) = r_x4e;
    r_y4e2 = R_y + rho_y4e2;  %y_record5(i) = r_y4e;
    r_z4e2 = R_z + rho_z4e2 + (3^0.5)/2*h;  %z_record5(i) = r_z4e;
    
    x_record6(i) = (r_x2e2+r_x4e2)/2;
    y_record6(i) = (r_y2e2+r_y4e2)/2;
    z_record6(i) = (r_z2e2+r_z4e2)/2;
    
    plotcircle1 = plot3(r_x,r_y,r_z,'b','linewidth',2);
    plotcircle2 = plot3(r2_x,r2_y,r2_z,'g','linewidth',2);
    plotcircle3 = plot3(r2_xe,r2_ye,r2_ze,'r','linewidth',2);
    plotcircle4 = plot3(r2_xe2,r2_ye2,r2_ze2,'m','linewidth',2);
    plotline0 = plot3([R_x; r_xvt], [R_y; r_yvt], [R_z; r_zvt],'b','linewidth',2);
    plotline1 = plot3([r_x2; r_x4], [r_y2; r_y4], [r_z2; r_z4],'g','linewidth',2); %!!
    plotline2 = plot3([r_x2e; r_x4e], [r_y2e; r_y4e], [r_z2e; r_z4e],'r','linewidth',2);
%     plotline3 = plot3([r_x2e2; r_x4e2], [r_y2e2; r_y4e2], [r_z2e2; r_z4e2],'m','linewidth',2);
    plotline3 = plot3([r_x2e2; (r_x2e2+r_x4e2)/2], [r_y2e2; (r_y2e2+r_y4e2)/2], [r_z2e2; (r_z2e2+r_z4e2)/2],'m','linewidth',2);
    pause(0.02);
    if (i < n)
        delete(plotcircle1);
        delete(plotcircle2);
        delete(plotcircle3);
        delete(plotcircle4);
        delete(plotline0);
        delete(plotline1);
        delete(plotline2);
        delete(plotline3);
    end
    
end
plot3(x_record1, y_record1, z_record1, 'b-.','linewidth',1);
plot3(x_record3, y_record3, z_record3, 'g-.','linewidth',1);
plot3(x_record4, y_record4, z_record4, 'r-.','linewidth',1);
% plot3(x_record5, y_record5, z_record5, 'm-.','linewidth',1);
% plot3(x_record6, y_record6, z_record6, 'm-.','linewidth',1);