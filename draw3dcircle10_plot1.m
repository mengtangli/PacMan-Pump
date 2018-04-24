% Plot the frame transformation relationship
% 
% Mengtang Li
% Apr 23

clear; clf; close all;

% color defination
% yellow = [0.87 0.49 0];
% green = [0 0.5 0];
figure('OuterPosition', [50 50 1000 1000]);
grid minor; grid on; hold on;
% plot frame IJK, positive and negative
plotaxisIp = plot3([0 2], [0 0], [0 0],'k','linewidth',2);
plotaxisJp = plot3([0 0], [0 2], [0 0],'k','linewidth',2);
plotaxisKp = plot3([0 0], [0 0], [0 2],'k','linewidth',2);
plotaxisIn = plot3([0 -2], [0 0], [0 0],'k-.','linewidth',2);
plotaxisJn = plot3([0 0], [0 -2], [0 0],'k-.','linewidth',2);
plotaxisKn = plot3([0 0], [0 0], [0 -2],'k-.','linewidth',2);
% specify the lengh of axis
xlim([-2.5 2.5]); ylim([-2.5 2.5]); zlim([-2.5 2.5]);
% legend('x-axis', 'y-axis', 'z-axis');
text(2.0,0.05,0.05,'I','FontSize', 16);
text(0.05,2.0,0.05,'J','FontSize', 16);
text(0.05,0.05,2.0,'K','FontSize', 16);
daspect([1 1 1]);
view(24,24);

alpha = pi/6;
beta = pi/6;
R = 1;
L = (sqrt(3)-1)*2*R;
h = L;
angle = 0:0.01:pi/6;
n = size(angle,2);

angle1 = 0:0.01:2*pi;
plot3(L/2*cos(angle1), L/2*sin(angle1), 0*angle1, 'r-.', 'linewidth',2)
plot3(R*cos(angle1), R*sin(angle1), 0*angle1+R, 'r', 'linewidth',2)
pause(2);

for i = 1:1:n
    wt = angle(i);    
    Ainv = [ (3^(1/2)*cos(2*wt))/2,  (3^(1/2)*sin(2*wt))/2,    1/2;
        -sin(2*wt),               cos(2*wt),              0;
        -cos(2*wt)/2,            -sin(2*wt)/2,         3^(1/2)/2;];    
    Binv = [ cos(-wt) sin(-wt) 0;
        -sin(-wt) cos(-wt) 0;
        0       0       1;];
    % Now the C is correct    
    C = Binv*Ainv;    
    % distance between origin and center of rotating circle
    R_x = -L/2*cos(2*wt);
    R_y = -L/2*sin(2*wt);
    R_z = 0;   
    % the rotating circle 
    theta = 0:0.02:2*pi;
    rho_x = 2*R*(cos(theta)*C(1,1)+sin(theta)*C(2,1));
    rho_y = 2*R*(cos(theta)*C(1,2)+sin(theta)*C(2,2));
    rho_z = 2*R*(cos(theta)*C(1,3)+sin(theta)*C(2,3));      
    % vertex of normal vector of radius R in moving frame 
    rho_xvt = 0+0+2*R*C(3,1);
    rho_yvt = 0+0+2*R*C(3,2);
    rho_zvt = 0+0+2*R*C(3,3);  
    % 1st point (0 deg) in moving frame 
    rho_x1 = 2*R*(cos(0)*C(1,1)+sin(0)*C(2,1));
    rho_y1 = 2*R*(cos(0)*C(1,2)+sin(0)*C(2,2));
    rho_z1 = 2*R*(cos(0)*C(1,3)+sin(0)*C(2,3));
    % 2nd point (90 deg) in moving frame 
    rho_x2 = 2*R*(cos(pi/2)*C(1,1)+sin(pi/2)*C(2,1));
    rho_y2 = 2*R*(cos(pi/2)*C(1,2)+sin(pi/2)*C(2,2));
    rho_z2 = 2*R*(cos(pi/2)*C(1,3)+sin(pi/2)*C(2,3));
    % 1st point (0 deg) in ijk frame 
    rho_i1 = 2*R*(cos(0)*Ainv(1,1)+sin(0)*Ainv(2,1));
    rho_j1 = 2*R*(cos(0)*Ainv(1,2)+sin(0)*Ainv(2,2));
    rho_k1 = 2*R*(cos(0)*Ainv(1,3)+sin(0)*Ainv(2,3));
    % 2nd point (90 deg) in ijk frame 
    rho_i2 = 2*R*(cos(pi/2)*Ainv(1,1)+sin(pi/2)*Ainv(2,1));
    rho_j2 = 2*R*(cos(pi/2)*Ainv(1,2)+sin(pi/2)*Ainv(2,2));
    rho_k2 = 2*R*(cos(pi/2)*Ainv(1,3)+sin(pi/2)*Ainv(2,3));   
    % center of moving circle in fixed frame 
    r_x = R_x + rho_x;
    r_y = R_y + rho_y;
    r_z = R_z + rho_z;
    % vertex of normal vector of radius R in fixed frame 
    r_xvt = R_x + rho_xvt; 
    r_yvt = R_y + rho_yvt; 
    r_zvt = R_z + rho_zvt;  
    % 1st point (0 deg) in fixed frame
    r_x1 = R_x + rho_x1; 
    r_y1 = R_y + rho_y1; 
    r_z1 = R_z + rho_z1;  
    % 2nd point (90 deg) in fixed frame 
    r_x2 = R_x + rho_x2;  
    r_y2 = R_y + rho_y2;  
    r_z2 = R_z + rho_z2;  
    % 1st point (0 deg) in ijk frame
    r_i1 = R_x + rho_i1; 
    r_j1 = R_y + rho_j1; 
    r_k1 = R_z + rho_k1;  
    % 2nd point (90 deg) in ijk frame 
    r_i2 = R_x + rho_i2;  
    r_j2 = R_y + rho_j2;  
    r_k2 = R_z + rho_k2; 
    % plot each frame
    plotcircle1 = plot3(r_x,r_y,r_z,'b','linewidth',2);
    plotline0 = plot3([R_x; r_xvt], [R_y; r_yvt], [R_z; r_zvt],'b','linewidth',2);
    plotline1 = plot3([R_x; 0], [R_y; 0], [R_z; 0],'b-.','linewidth',2);
    plotline2 = plot3([R_x; r_x1], [R_y; r_y1], [R_z; r_z1],'m-.','linewidth',2);
    plotline3 = plot3([R_x; r_x2], [R_y; r_y2], [R_z; r_z2],'m-.','linewidth',2);
    plotline4 = plot3([R_x; r_i1], [R_y; r_j1], [R_z; r_k1],'b-.','linewidth',2);
    plotline5 = plot3([R_x; r_i2], [R_y; r_j2], [R_z; r_k2],'b-.','linewidth',2);
    pause(0.02);
    % and delete them after some pause
    if (i < n)
        delete(plotcircle1);
        delete(plotline0);
        delete(plotline1);
        delete(plotline2);
        delete(plotline3);
        delete(plotline4);
        delete(plotline5);       
    end
    
end
text(r_i1+0.1,r_j1+0.1,r_k1+0.1,'i','FontSize', 16);
text(r_i2+0.1,r_j2+0.1,r_k2+0.1,'j','FontSize', 16);
text(r_xvt+0.1,r_yvt+0.1,r_zvt+0.1,'k(z)','FontSize', 16);
text(r_x1+0.1,r_y1+0.1,r_z1+0.1,'x','FontSize', 16);
text(r_x2+0.1,r_y2+0.1,r_z2+0.1,'y','FontSize', 16);
% text(r_xvt+0.05,r_yvt+0.05,r_zvt+0.05,'z');