% ---- test 1 ---- %
phi = 0:0.01:2*pi;
x = 10+8*cos(phi)+3*sqrt(2)*sin(phi);
y = 20-6*cos(phi)+4*sqrt(2)*sin(phi);
z = 30-5*sqrt(2)*sin(phi);

figure(1);
plot3(x,y,z);
grid on; grid minor;
clear;

% ---- test 2 ---- %
phi = 0:0.01:2*pi;
alpha = pi/6;
R = 1;
x = R*cos(alpha)*sin(phi);
y = -R*cos(phi);
z = -R*sin(alpha)*sin(phi);

figure(2);
plot3(x,y,z,'k','linewidth',2);
grid on; grid minor; hold on;

xaxis_x = 0:0.01:2;
xaxis_y = zeros(1,size(xaxis_x,2));
xaxis_z = zeros(1,size(xaxis_x,2));

yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.01:2;
yaxis_z = zeros(1,size(xaxis_x,2));

zaxis_x = zeros(1,size(xaxis_x,2));
zaxis_y = zeros(1,size(xaxis_x,2));
zaxis_z = 0:0.01:2;

plot3(xaxis_x,xaxis_y,xaxis_z,'b','linewidth',2);
plot3(yaxis_x,yaxis_y,yaxis_z,'r','linewidth',2);
plot3(zaxis_x,zaxis_y,zaxis_z,'g','linewidth',2);
clear;

% ---- test 3 ---- %
phi = 0:0.01:2*pi;
beta = pi/6;
R = 1;
x = -R*sin(beta)*sin(2*phi);
y = -R*sin(beta)*cos(2*phi);
z = zeros(1,size(x,2));

figure(3);
plot3(x,y,z,'k','linewidth',2);
grid on; grid minor; hold on;

xaxis_x = 0:0.01:2;
xaxis_y = zeros(1,size(xaxis_x,2));
xaxis_z = zeros(1,size(xaxis_x,2));

yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.01:2;
yaxis_z = zeros(1,size(xaxis_x,2));

zaxis_x = zeros(1,size(xaxis_x,2));
zaxis_y = zeros(1,size(xaxis_x,2));
zaxis_z = 0:0.01:2;

plot3(xaxis_x,xaxis_y,xaxis_z,'b','linewidth',2);
plot3(yaxis_x,yaxis_y,yaxis_z,'r','linewidth',2);
plot3(zaxis_x,zaxis_y,zaxis_z,'g','linewidth',2);
clear;