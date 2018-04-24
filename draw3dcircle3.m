% Try to calculate the distance
% Mengtang Li
% Mar 25

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

xlim([-1.2 1.2]); ylim([-1.2 1.2]); zlim([-1.2 1.2]);

plot3(xaxis_x,xaxis_y,xaxis_z,'k','linewidth',2);
grid on; grid minor; hold on;
plot3(yaxis_x,yaxis_y,yaxis_z,'k','linewidth',2);
plot3(zaxis_x,zaxis_y,zaxis_z,'k','linewidth',2);
% legend('x-axis', 'y-axis', 'z-axis');
text(1.0,0.05,0.05,'x');
text(0.05,1.0,0.05,'y');
text(0.05,0.05,1.0,'z');
daspect([1 1 1]);

alpha = pi/6;
beta = pi/6;
R = 1;
wt = pi/4;

% lower
x1 = [0; -R*cos(alpha)*sin(wt)];
y1 = [0; -R*cos(wt)];
z1 = [0; R*sin(alpha)*sin(wt)];
plot3(x1, y1, z1,'b','linewidth',2);

% upper
theta = 0:0.05:2*pi;
x2 = R*( -sin(2*wt)*cos(theta) + cos(beta)*cos(2*wt)*sin(theta) );
y2 = R*( -cos(2*wt)*cos(theta) - cos(beta)*sin(2*wt)*sin(theta) );
z2 = R*(  sin(beta)*sin(theta) );
plot3(x2, y2, z2,'r','linewidth',2);
% upper normal
x2_n = [0; R*( -sin(beta)*cos(2*wt) )];
y2_n = [0; R*(  sin(beta)*sin(2*wt) )];
z2_n = [0; R*(  cos(beta) )];
plot3(x2_n, y2_n, z2_n,'r','linewidth',2);

% possible upper closest point 1
x2_c = [0; R*( -sin(2*wt)*cos(wt) + cos(beta)*cos(2*wt)*sin(wt) )];
y2_c = [0; R*( -cos(2*wt)*cos(wt) - cos(beta)*sin(2*wt)*sin(wt) )];
z2_c = [0; R*( sin(beta)*sin(wt) )];
plot3(x2_c, y2_c, z2_c,'r-.','linewidth',2);

