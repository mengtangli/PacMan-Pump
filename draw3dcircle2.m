% By changing (guessing) the signs of x1 y1 z1 x2 y2 z2 to make the
% rotation direction consistent with each other.
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

xlim([-1 1.5]); ylim([-1 1.5]); zlim([-1 1.5]);

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

% lower
for wt = 0:pi/8:3*pi/4;
    x1 = R*cos(alpha)*sin(wt);
    y1 = -R*cos(wt);
    z1 = -R*sin(alpha)*sin(wt);
    x1 = [0; x1];
    y1 = [0; y1];
    z1 = [0; z1];
    plot3(x1, y1, z1,'b','linewidth',2);
end

% upper
wt = 0;
for theta = 0:pi/8:3*pi/4;
    x2 = R*( -sin(2*wt)*cos(theta) + cos(beta)*cos(2*wt)*sin(theta) );
    y2 = R*( -cos(2*wt)*cos(theta) - cos(beta)*sin(2*wt)*sin(theta) );
    z2 = R*( +sin(beta)*sin(theta) );
    x2 = [0; x2];
    y2 = [0; y2];
    z2 = [0; z2];
    plot3(x2, y2, z2,'r','linewidth',2);
end