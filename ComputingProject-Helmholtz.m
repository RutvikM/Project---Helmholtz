clear all,clc
Ax = -pi; % Left most point
Bx = pi;
Ay = -pi;
By = pi;
pt = 10; % The number of points, NOT the number of discretizations.
dx = (Bx-Ax)/pt;
dy = (By-Ay)/pt;
x = zeros(1,pt);
y = zeros(1,pt);
x(1,1) = Ax;
y(1,1) = Ay;
for a = 2:pt
    x(1,a) = Ax+(a*dx); % Grid being made in the for loop.
    y(1,a) = Ay+(a*dx);
end
