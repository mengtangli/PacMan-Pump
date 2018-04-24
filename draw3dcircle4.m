% Dont know what I am going to draw
% http://sycomoreen.free.fr/docs_multimedia/Double_Spherical_Wankel16.html
% Run one section from the following three
% Mengtang Li
% Mar 28
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


% % ------------- Plot one Circle and points ------------- %
% wt = pi;
% for theta = 0:pi/4:pi;
%     x = [0; R*( -sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) )];
%     y = [0; R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) )];
%     z = [0; R*( -sin(alpha)*sin(theta) )];
%     plot3(x, y, z,'r','linewidth',2);
% end
% theta = 0:0.01:2*pi;
% x = R*( -sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) );
% y = R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) );
% z = R*( -sin(alpha)*sin(theta) );
% 
% plot3(x,y,z,'r-','linewidth',2);
% clear x y z;


% % ------------- Plot Rotating Cirlce ------------- %
% theta = 0:0.01:2*pi;
% % ---- wt = 0 ---- %
% wt = 0;
% x = R*( -sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) );
% y = R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) );
% z = R*( -sin(alpha)*sin(theta) );
% plot3(x,y,z,'k-','linewidth',2);
% clear x y z;
% % ---- wt = pi/4 ---- %
% wt = pi/4;
% x = R*( -sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) );
% y = R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) );
% z = R*( -sin(alpha)*sin(theta) );
% plot3(x,y,z,'r-','linewidth',2);
% clear x y z;
% % ---- wt = pi/2 ---- %
% wt = pi/2;
% x = R*( -sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) );
% y = R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) );
% z = R*( -sin(alpha)*sin(theta) );
% plot3(x,y,z,'r-.','linewidth',2);
% clear x y z;
% % ---- wt = 3pi/2 ---- %
% wt = 3*pi/2;
% x = R*( sin(wt)*cos(theta) + cos(alpha)*cos(wt)*sin(theta) );
% y = R*( -cos(wt)*cos(theta) - cos(alpha)*sin(wt)*sin(theta) );
% z = R*( -sin(alpha)*sin(theta) );
% plot3(x,y,z,'r:','linewidth',2);
% clear x y z;


% ------------- Plot Trajectory ------------- %
wt = 0:0.01:2*pi;
x1 = R*( -sin(wt).*cos(wt) + cos(alpha)*cos(wt).*sin(wt) );
y1 = R*( -cos(wt).*cos(wt) - cos(alpha)*sin(wt).*sin(wt) );
z1 = R*( sin(wt)*0 - sin(alpha)*sin(wt) );
plot3(x1,y1,z1,'r','linewidth',2);

x2 = R*( -sin(wt).*cos(wt+pi) + cos(alpha)*cos(wt).*sin(wt+pi) );
y2 = R*( -cos(wt).*cos(wt+pi) - cos(alpha)*sin(wt).*sin(wt+pi) );
z2 = R*( sin(wt)*0 - sin(alpha)*sin(wt+pi) );
plot3(x2,y2,z2,'r','linewidth',2);
