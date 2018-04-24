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

xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);

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
L = 2;
angle = 0:0.05:2*pi;
n = size(angle,2);
x_record = zeros(1,n);
y_record = zeros(1,n);
z_record = zeros(1,n);

wt = 0
A = [ (3^(1/2)*cos(wt))/(2*(cos(wt)^2 + sin(wt)^2)), (3^(1/2)*sin(wt))/(2*(cos(wt)^2 + sin(wt)^2)),    -1/2;
    -sin(wt)/(cos(wt)^2 + sin(wt)^2),               cos(wt)/(cos(wt)^2 + sin(wt)^2),                    0;
    cos(wt)/(2*(cos(wt)^2 + sin(wt)^2)),           sin(wt)/(2*(cos(wt)^2 + sin(wt)^2)),       3^(1/2)/2;];

x = [0 R*0*A(1,1)-R*1*A(2,1)+L*A(3,1)];
y = [0 R*0*A(1,2)-R*1*A(2,2)+L*A(3,2)];
z = [0 R*0*A(1,3)-R*1*A(2,3)+L*A(3,3)-sqrt(3)/2*L];

plot3(x,y,z,'r','linewidth',2);

% for i = 1:1:n
%     wt = angle(i);
%     A = [ (3^(1/2)*cos(wt))/(2*(cos(wt)^2 + sin(wt)^2)), (3^(1/2)*sin(wt))/(2*(cos(wt)^2 + sin(wt)^2)),    -1/2;
%         -sin(wt)/(cos(wt)^2 + sin(wt)^2),               cos(wt)/(cos(wt)^2 + sin(wt)^2),                    0;
%         cos(wt)/(2*(cos(wt)^2 + sin(wt)^2)),           sin(wt)/(2*(cos(wt)^2 + sin(wt)^2)),       3^(1/2)/2;];
%     
%     %     x = R*sin(wt)*A(1,1)-R*cos(wt)*A(2,1)+L*A(3,1);
%     %     y = R*sin(wt)*A(1,2)-R*cos(wt)*A(2,2)+L*A(3,2);
%     %     z = R*sin(wt)*A(1,3)-R*cos(wt)*A(2,3)+L*A(3,3);%-sqrt(3)/2*L;
%     x = R*0*A(1,1)-R*1*A(2,1)+L*A(3,1);
%     y = R*0*A(1,2)-R*1*A(2,2)+L*A(3,2);
%     z = R*0*A(1,3)-R*1*A(2,3)+L*A(3,3)-sqrt(3)/2*L;
%     
%     x_record(i) = x;
%     y_record(i) = y;
%     z_record(i) = z;
% end
% 
% plot3(x_record,y_record,z_record,'r','linewidth',2);