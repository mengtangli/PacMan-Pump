% Run this after runing draw3dcircle10_test4.m
% Mengtang Li
% mar 13

figure(1);
hold on;


angle = 0:0.02:4*pi;
n = size(angle,2);
x_re1 = zeros(1,n);
y_re1 = zeros(1,n);
z_re1 = zeros(1,n);
x_re2 = zeros(1,n);
y_re2 = zeros(1,n);
z_re2 = zeros(1,n);
R = 1;
L = (sqrt(3)-1)*2*R;
h = L;

% for i = 1:1:n
%     t = angle(i);
%     R_x = -L/2*cos(2*t);
%     R_y = -L/2*sin(2*t);
%     R_z = 0;
%     rho_x1 = 3*cos(pi+0.1002)*(sqrt(3)/2*sin(t)*cos(2*t)-cos(t)*sin(2*t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(-1/2*cos(2*t));
%     rho_y1 = 3*cos(pi+0.1002)*(sqrt(3)/2*sin(t)*sin(2*t)+cos(t)*cos(2*t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(-1/2*sin(2*t));
%     rho_z1 = 3*cos(pi+0.1002)*(1/2*sin(t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(sqrt(3)/2);
%     rho_x2 = 3*cos(pi-0.1002)*(sqrt(3)/2*sin(t)*cos(2*t)-cos(t)*sin(2*t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(-1/2*cos(2*t));
%     rho_y2 = 3*cos(pi-0.1002)*(sqrt(3)/2*sin(t)*sin(2*t)+cos(t)*cos(2*t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(-1/2*sin(2*t));
%     rho_z2 = 3*cos(pi-0.1002)*(1/2*sin(t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(sqrt(3)/2);
%     r_x1 = R_x + rho_x1;
%     r_y1 = R_y + rho_y1;
%     r_z1 = R_z + rho_z1 + (3^0.5)/2*h;
%     r_x2 = R_x + rho_x2;
%     r_y2 = R_y + rho_y2;
%     r_z2 = R_z + rho_z2 + (3^0.5)/2*h;
%     x_re1(i) = r_x1;
%     y_re1(i) = r_y1;
%     z_re1(i) = r_z1;
%     x_re2(i) = r_x2;
%     y_re2(i) = r_y2;
%     z_re2(i) = r_z2;
% end

for i = 1:1:n
    t = angle(i);
    R_x = -L/2*cos(2*t);
    R_y = -L/2*sin(2*t);
    R_z = 0;  
    rho_x1 = 3*cos(pi+0.1002)*(sqrt(3)/2*sin(t)*cos(2*t)-cos(t)*sin(2*t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(-1/2*cos(2*t)) -(sqrt(3)-1)*cos(2*t);
    rho_y1 = 3*cos(pi+0.1002)*(sqrt(3)/2*sin(t)*sin(2*t)+cos(t)*cos(2*t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(-1/2*sin(2*t)) -(sqrt(3)-1)*sin(2*t);
    rho_z1 = 3*cos(pi+0.1002)*(1/2*sin(t)) + (3*sin(pi+0.1002)-2*(sqrt(3)-1))*(sqrt(3)/2) + (3-sqrt(3));
    rho_x2 = 3*cos(pi-0.1002)*(sqrt(3)/2*sin(t)*cos(2*t)-cos(t)*sin(2*t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(-1/2*cos(2*t));
    rho_y2 = 3*cos(pi-0.1002)*(sqrt(3)/2*sin(t)*sin(2*t)+cos(t)*cos(2*t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(-1/2*sin(2*t));
    rho_z2 = 3*cos(pi-0.1002)*(1/2*sin(t)) + (3*sin(pi-0.1002)-2*(sqrt(3)-1))*(sqrt(3)/2);
    r_x1 = 0 + rho_x1;
    r_y1 = 0 + rho_y1;
    r_z1 = 0 + rho_z1;
    r_x2 = R_x + rho_x2;
    r_y2 = R_y + rho_y2;
    r_z2 = R_z + rho_z2 + (3^0.5)/2*h;
    x_re1(i) = r_x1;
    y_re1(i) = r_y1;
    z_re1(i) = r_z1;
    x_re2(i) = r_x2;
    y_re2(i) = r_y2;
    z_re2(i) = r_z2;
end

plot3(x_re1, y_re1, z_re1, 'b-.','linewidth',2);
plot3(x_re2, y_re2, z_re2, 'b-.','linewidth',2);
hold off;