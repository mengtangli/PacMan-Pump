% plot out the trajectory of 3D gear set 2 vertex (if 2, I guess it doesnt matter 
% what the gear ratio is?), no L
% Mengtang Li
% Apr 10
clear; clf; close all;

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
L = (2*sqrt(3)-3)*R;
h = L;
angle = 0:0.02:4*pi;
n = size(angle,2);
x_record1 = zeros(1,n);
y_record1 = zeros(1,n);
z_record1 = zeros(1,n);
x_record2 = zeros(1,n);
y_record2 = zeros(1,n);
z_record2 = zeros(1,n);

plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle+sqrt(7)/2*R, 'r', 'linewidth',2)
plot3(1.5*R*cos(angle), 1.5*R*sin(angle), 0*angle-sqrt(7)/2*R, 'r', 'linewidth',2)
pause(2);

for i = 1:1:n
    wt = angle(i);
    
    x_record1(i) = 3*(sqrt(3)/2*cos(wt)*cos(2*wt)+sin(wt)*sin(2*wt))+h*1/2*cos(2*wt)-L/2*cos(2*wt);
    y_record1(i) = 3*(sqrt(3)/2*cos(wt)*sin(2*wt)-sin(wt)*cos(2*wt))+h*1/2*sin(2*wt)-L/2*sin(2*wt);
    z_record1(i) = 3*(1/2*cos(wt))-sqrt(3)/2*h+sqrt(3)/2*L;
    
    x_record2(i) = 3*(sqrt(3)/2*sin(wt)*cos(2*wt)-cos(wt)*sin(2*wt))+h*1/2*cos(2*wt)-L/2*cos(2*wt);
    y_record2(i) = 3*(sqrt(3)/2*sin(wt)*sin(2*wt)+cos(wt)*cos(2*wt))+h*1/2*sin(2*wt)-L/2*sin(2*wt);
    z_record2(i) = 3*(1/2*cos(wt))-sqrt(3)/2*h+sqrt(3)/2*L;   
end
plot3(x_record1, y_record1, z_record1, 'r-.','linewidth',1);
plot3(x_record2, y_record2, z_record2, 'b-.','linewidth',1);
hold off;