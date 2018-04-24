% Try to use vector method to draw the circle
% Mengtang Li
% Mar 30
clear; clf; close all;

figure(1);
grid minor; grid on; hold on;
% plot axes
xaxis_x = 0:0.1:1;
xaxis_y = zeros(1,size(xaxis_x,2));
xaxis_z = zeros(1,size(xaxis_x,2));

yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.1:1;
yaxis_z = zeros(1,size(xaxis_x,2));

zaxis_x = zeros(1,size(xaxis_x,2));
zaxis_y = zeros(1,size(xaxis_x,2));
zaxis_z = 0:0.1:1;

% xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

plot3(xaxis_x,xaxis_y,xaxis_z,'k','linewidth',2);
grid on; grid minor; hold on;
plot3(yaxis_x,yaxis_y,yaxis_z,'k','linewidth',2);
plot3(zaxis_x,zaxis_y,zaxis_z,'k','linewidth',2);
% legend('x-axis', 'y-axis', 'z-axis');
text(1.0,0.05,0.05,'x');
text(0.05,1.0,0.05,'y');
text(0.05,0.05,1.0,'z');
daspect([1 1 1]);
view(24,24);

alpha = pi/6;
beta = pi/6;
R = 1;
L = 2;
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

for i = 1:1:n
    wt = angle(i);
%     A = [ (3^(1/2)*cos(wt))/(2*(cos(wt)^2 + sin(wt)^2)), (3^(1/2)*sin(wt))/(2*(cos(wt)^2 + sin(wt)^2)),    -1/2;
%         -sin(wt)/(cos(wt)^2 + sin(wt)^2),               cos(wt)/(cos(wt)^2 + sin(wt)^2),                    0;
%         cos(wt)/(2*(cos(wt)^2 + sin(wt)^2)),           sin(wt)/(2*(cos(wt)^2 + sin(wt)^2)),       3^(1/2)/2;];
    
    A = [ (3^(1/2)*cos(wt))/2, (3^(1/2)*sin(wt))/2,    -1/2;
        -sin(wt),               cos(wt),                    0;
        cos(wt)/2,           sin(wt)/2,       3^(1/2)/2;];
    
    rho_x1 = R*sin(wt)*A(1,1)-R*cos(wt)*A(2,1)+L*A(3,1);
    rho_y1 = R*sin(wt)*A(1,2)-R*cos(wt)*A(2,2)+L*A(3,2);
    rho_z1 = R*sin(wt)*A(1,3)-R*cos(wt)*A(2,3)+L*A(3,3);
    rho_x2 = R*sin(wt+pi)*A(1,1)-R*cos(wt+pi)*A(2,1)+L*A(3,1);
    rho_y2 = R*sin(wt+pi)*A(1,2)-R*cos(wt+pi)*A(2,2)+L*A(3,2);
    rho_z2 = R*sin(wt+pi)*A(1,3)-R*cos(wt+pi)*A(2,3)+L*A(3,3);
    R_x = -L/2*cos(wt);
    R_y = -L/2*sin(wt);
    R_z = -sqrt(3)/2*L;
    
    x_lower = -R*cos(alpha)*sin(wt);
    y_lower = -R*cos(wt);
    z_lower = R*sin(alpha)*sin(wt);
    
    x_record1(i) = rho_x1+R_x;
    y_record1(i) = rho_y1+R_y;
    z_record1(i) = rho_z1+R_z;
    x_record2(i) = rho_x2+R_x;
    y_record2(i) = rho_y2+R_y;
    z_record2(i) = rho_z2+R_z;
    x_record3(i) = x_lower;
    y_record3(i) = y_lower;
    z_record3(i) = z_lower;
end

plot3(x_record1,y_record1,z_record1,'r-.','linewidth',2);
plot3(x_record2,y_record2,z_record2,'b-.','linewidth',2);
plot3(x_record3,y_record3,z_record3,'g-.','linewidth',2);

for i = 1:1:n
    hh1 = plot3([0; x_record1(i)], [0; y_record1(i)], [0; z_record1(i)], 'k', 'linewidth',2);
    hh2 = plot3([0; x_record2(i)], [0; y_record2(i)], [0; z_record2(i)], 'k', 'linewidth',2);
    pause(0.02);
    delete(hh1);
    delete(hh2);
end