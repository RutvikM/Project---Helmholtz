% SOR
clear all,clc
Ax = -pi;
Bx = pi;
Ay = -pi;
By = pi;
w = 1.4;
N = 4; % Number of points; N+1 discretizations % ONLY WORKS FOR UP TO 4.
h = (Bx-Ax)/(N+2); % Since the discretization for x and y are equal, h is used.
x = zeros(1,(N+2));
y = zeros(1,(N+2));
x(1,1) = Ax;
y(1,1) = Ay;
for a = 2:(N+2)
    x(1,a) = Ax+(a*h); % Grid being made in the for loop.
    y(1,a) = Ay+(a*h);
end
counter = 1;
epsilon = 0.1;
Error = 10;
U = ones(N+2,N+2); % This is the initial guess matrix.
C = ones(N+2,N+2); % This is the matrix that will copy the previous time frame matrix 
for i = 1:(N+2)
        C(1,i) = cos(pi*(x(i)-Ax))*cosh(Bx-x(i));
        C((N+2),i) = (x(i)-Ax)^2*sin(pi*(x(i)-Ax)/(2*(Bx-Ax)));
end
while Error>epsilon
    for i = 2:(N+1)
       for j = 2:(N+1)
             F(i,j) = cos(pi/2.*(2.*((x(j)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i)-Ay)/(By-Ay)));
              C(i,j) = 0.25*w*(C(i-1,j)+C(i,j-1)+U(i,j+1)+U(i+1,j)-(F(i,j)*h^2))+((1-w)*U(i,j));
       end
         F(i,1) = cos(pi/2.*(2.*((x(1)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i))-Ay)/(By-Ay)); 
         C(i,1) = 0.25*w*(C(i-1,1)+C(i,2)+U(i,2)+U(i+1,1)-(F(i,1)*h^2))+((1-w)*U(i,j));
         F(i,N+2) = cos(pi/2.*(2.*((x(N+2)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i))-Ay)/(By-Ay)); 
         C(i,(N+2)) = 0.25*w*(C(i-1,N+2)+C(i,N+1)+U(i,N+1)+U(i+1,N+2)-(h^2*F(i,N+2)))+((1-w)*U(i,j));
    end
Err = max(abs(C-U));
Error = max(Err);
U = C;
counter = counter+1;
end
surf(x,y,C)
title(['Solution of Poisson Equation Using SOR - Iteration ',num2str(counter)]);
xlabel('X','Fontsize',14)
ylabel('Y','Fontsize',14)
zlabel('U','Fontsize',14)