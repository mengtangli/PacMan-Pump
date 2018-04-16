% Plot moving 2D gear set
% Mengtang Li
% Mar 30
clear; clf; close all;

figure('OuterPosition', [50 50 1000 1000]);
grid minor; grid on; hold on;
% plot axes
xaxis_x = 0:0.1:2;
xaxis_y = zeros(1,size(xaxis_x,2));

yaxis_x = zeros(1,size(xaxis_x,2));
yaxis_y = 0:0.1:2;

xlim([-3 3]); ylim([-3 3]);

plot(xaxis_x,xaxis_y,'k','linewidth',2);
grid on; grid minor; hold on;
plot(yaxis_x,yaxis_y,'k','linewidth',2);
% legend('x-axis', 'y-axis', 'z-axis');
text(2.0,0.05,'x');
text(0.05,2.0,'y');
daspect([1 1 1]);

alpha = pi/6;
beta = pi/6;
R = 1;
L = 2;
angle = 0:0.02:2*pi;
n = size(angle,2);
x_record = zeros(1,n);
y_record = zeros(1,n);

circle1 = viscircles([0 0], R);

for i = 1:1:n
    wt = angle(i);
    %     A = [ cos(wt)/(cos(wt)^2 + sin(wt)^2),  sin(wt)/(cos(wt)^2 + sin(wt)^2);
    %         -sin(wt)/(cos(wt)^2 + sin(wt)^2),  cos(wt)/(cos(wt)^2 + sin(wt)^2);];
    
    A = [ cos(wt),  sin(wt);
        -sin(wt),  cos(wt);]; % bigger gear also rotates
    
    %     A = [1 0; 0 1]; % bigger gear doesnt rotate
    
    R_x = R*cos(2*wt);
    R_y = R*sin(2*wt);
    
    theta = 0:0.02:2*pi;
    rho_x = 2*R*(-sin(theta)*A(1,1)+cos(theta)*A(2,1));
    rho_y = 2*R*(-sin(theta)*A(1,2)+cos(theta)*A(2,2));
    
    rho_x1 = 2*R*(-sin(0)*A(1,1)+cos(0)*A(2,1));
    rho_y1 = 2*R*(-sin(0)*A(1,2)+cos(0)*A(2,2));
    rho_x2 = 2*R*(-sin(pi/2)*A(1,1)+cos(pi/2)*A(2,1));
    rho_y2 = 2*R*(-sin(pi/2)*A(1,2)+cos(pi/2)*A(2,2));
    rho_x3 = 2*R*(-sin(pi)*A(1,1)+cos(pi)*A(2,1));
    rho_y3 = 2*R*(-sin(pi)*A(1,2)+cos(pi)*A(2,2));
    rho_x4 = 2*R*(-sin(3*pi/2)*A(1,1)+cos(3*pi/2)*A(2,1));
    rho_y4 = 2*R*(-sin(3*pi/2)*A(1,2)+cos(3*pi/2)*A(2,2));
    
    r_x = R_x + rho_x;
    r_y = R_y + rho_y;
    
    r_x1 = R_x + rho_x1;
    r_y1 = R_y + rho_y1;
    r_x2 = R_x + rho_x2;
    r_y2 = R_y + rho_y2;
    r_x3 = R_x + rho_x3;
    r_y3 = R_y + rho_y3;
    r_x4 = R_x + rho_x4;
    r_y4 = R_y + rho_y4;
    
    plotcircle = plot(r_x,r_y,'b','linewidth',2);
    plotline1 = plot([r_x1; r_x3], [r_y1; r_y3], 'g','linewidth',2);
    plotline2 = plot([r_x2; r_x4], [r_y2; r_y4], 'y','linewidth',2);
    
    pause(0.1);
    if (i < n)
        delete(plotcircle);
        delete(plotline1);
        delete(plotline2);
    end
    
end
