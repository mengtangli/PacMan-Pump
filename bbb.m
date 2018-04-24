% check whether 10_test5 is correct
% Apr 14
% Mengtang


% alpha = pi/6;
% beta = pi/6;
% R = 1;
% L = (sqrt(3)-1)*2*R;
% h = L;

alpha = 9*pi/180;
beta = 18*pi/180;
R = 1;
%L = R*(2*sin(beta)-tan(beta))/(tan(beta)*sin(beta));
L = R*(2*cos(beta)-1)/sin(beta);
h = L;

angle = 0:0.01:2*pi;
n = size(angle,2);
x_record0 = zeros(1,n);
y_record0 = zeros(1,n);
z_record0 = zeros(1,n);

for i = 1:1:n
    wt = angle(i);

%     C11 = sqrt(3)/2*cos(wt)*cos(2*wt)+sin(wt)*sin(2*wt);
%     C12 = sqrt(3)/2*cos(wt)*sin(2*wt)-sin(wt)*cos(2*wt);
%     C13 = 1/2*cos(wt);
%     C21 = sqrt(3)/2*sin(wt)*cos(2*wt)-cos(wt)*sin(2*wt);
%     C22 = sqrt(3)/2*sin(wt)*sin(2*wt)+cos(wt)*cos(2*wt);
%     C23 = 1/2*sin(wt);
%     C31 = -1/2*cos(2*wt);
%     C32 = -1/2*sin(2*wt);
%     C33 = sqrt(3)/2;
%     R_x = -L/2*cos(2*wt)+0.3*R*0.5;
%     R_y = -L/2*sin(2*wt)+0;
%     R_z = 0+(3^0.5)/2*h-0.3*R*0.5*sqrt(3);

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
    R_x = -L*sin(beta)*cos(2*wt);
    R_y = -L*sin(beta)*sin(2*wt);
    R_z = cos(beta)*h;

    rho_x = 3*R*cos(pi+asin(0.1))*C21 + (3*R*sin(pi+asin(0.1))-h)*C31;
    rho_y = 3*R*cos(pi+asin(0.1))*C22 + (3*R*sin(pi+asin(0.1))-h)*C32;
    rho_z = 3*R*cos(pi+asin(0.1))*C23 + (3*R*sin(pi+asin(0.1))-h)*C33;

    x_record0(i) = R_x+rho_x;
    y_record0(i) = R_y+rho_y;
    z_record0(i) = R_z+rho_z;

end

% for i = 1:1:n
%     wt = angle(i);
%     % wrong?
% %     x_record0(i) = -(sqrt(3)-1)*cos(2*wt)+2.5 + 3*cos(pi+asin(0.1))*(sqrt(3)/2*sin(wt)*cos(2*wt)-cos(wt)*sin(2*wt)) + (3*sin(pi+asin(0.1))-(2*sqrt(3)-7))*(-1/2*cos(2*wt));
% %     y_record0(i) = -(sqrt(3)-1)*sin(2*wt) + 3*cos(pi+asin(0.1))*(sqrt(3)/2*sin(wt)*sin(2*wt)+cos(wt)*cos(2*wt)) + (3*sin(pi+asin(0.1))-(2*sqrt(3)-7))*(-1/2*sin(2*wt));
% %     z_record0(i) = (3-sqrt(3))-2.5*sqrt(3) + 3*cos(pi+asin(0.1))*(1/2*sin(wt)) + (3*sin(pi+asin(0.1))-(2*sqrt(3)-7))*(sqrt(3)/2);
%     % correct
%     x_record0(i) = -(sqrt(3)-1)*cos(2*wt)+0.3*R*0.5 + 3*R*cos(pi+asin(0.1))*(sqrt(3)/2*sin(wt)*cos(2*wt)-cos(wt)*sin(2*wt)) + (3*R*sin(pi+asin(0.1))-(h-3*R*0.1))*(-1/2*cos(2*wt));
%     y_record0(i) = -(sqrt(3)-1)*sin(2*wt) + 3*R*cos(pi+asin(0.1))*(sqrt(3)/2*sin(wt)*sin(2*wt)+cos(wt)*cos(2*wt)) + (3*R*sin(pi+asin(0.1))-(h-3*R*0.1))*(-1/2*sin(2*wt));
%     z_record0(i) = (3-sqrt(3))-0.3*R*0.5*sqrt(3) + 3*R*cos(pi+asin(0.1))*(1/2*sin(wt)) + (3*R*sin(pi+asin(0.1))-(h-3*R*0.1))*(sqrt(3)/2);
%     
% end
figure(1); hold on;
plot3(x_record0, y_record0, z_record0, 'b:','linewidth',2);
hold off;
