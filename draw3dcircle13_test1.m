
clear;clf;


R = 1;
angle = 0:0.02:4*pi;
n = size(angle,2);
x_record1 = zeros(1,n);
y_record1 = zeros(1,n);
z_record1 = zeros(1,n);

figure(1);
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
    
    rho1 = C'*[cos(0); sin(0); 0]*2*R;
    R1 = [0; 0; 0];
    r1 = R1+rho1;
    x_record1(i) = r1(1);
    y_record1(i) = r1(2);
    z_record1(i) = r1(3);
end

plot3(x_record1, y_record1, z_record1, 'r-.','linewidth',1);
hold off;