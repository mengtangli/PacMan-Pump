% Calculate the volume: find the angle between normal vectors.
% Continued after 14_3.m
% Mengtang Li
% Apr 20

clear;clf;

alpha = 9*pi/180;
beta = 18*pi/180;
R = 1;
%L = R*(2*sin(beta)-tan(beta))/(tan(beta)*sin(beta));
L = R*(2*cos(beta)-1)/sin(beta);
h = L;
r1 = 50;
r2 = 30;

angle = 0:0.01:2*pi;
n = size(angle,2);
angle_record = zeros(1,n);

for i = 1:1:n
    wt = angle(i);
    
    % ------- General angle: beta
    C11 = cos(beta)*cos(wt)*cos(2*wt)+sin(wt)*sin(2*wt);
    C12 = cos(beta)*cos(wt)*sin(2*wt)-sin(wt)*cos(2*wt);
    C13 = sin(beta)*cos(wt);
    C21 = cos(beta)*sin(wt)*cos(2*wt)-cos(wt)*sin(2*wt);
    C22 = cos(beta)*sin(wt)*sin(2*wt)+cos(wt)*cos(2*wt);
    C23 = sin(beta)*sin(wt);
    C31 = -sin(beta)*cos(2*wt);
    C32 = -sin(beta)*sin(2*wt);
    C33 = cos(beta);
    
    n1 = [3*R*sin(beta) 0 3*R*cos(beta)];
    n2 = [-h*sin(beta)*cos(2*wt)+3*R*C11-h*C31;
        -h*sin(beta)*sin(2*wt)+3*R*C12-h*C32;
        h*cos(beta)+3*R*C13-h*C33;];
    angle_record(i) = acos(n1*n2/(norm(n1)*norm(n2)));
end
figure(1);
plot(angle*180/pi, angle_record*180/pi, 'b', 'linewidth', 2);
grid minor; grid on;
xlabel('rotated angle, deg');
ylabel('volume angle, deg');
hold on;
plot(angle*180/pi, (126-54)/2*sin(angle-pi/2)+(126+54)/2, 'r-.', 'linewidth', 2);
hold off;
legend('calculation', 'pure sine');

figure(2);
Vmin = 4/3*pi*45/360*(r1^3-r2^3)/1000;
Vmax = 4/3*pi*117/360*(r1^3-r2^3)/1000;
plot(angle*180/pi, (Vmax-Vmin)/2*sin(angle-pi/2)+(Vmax+Vmin)/2, 'b', 'linewidth', 2);
grid minor; grid on;
xlabel('rotated angle, deg');
ylabel('volume, mL');
