%% USING GAUSS SEIDEL METHOD
% PART A
clear all,clc
Ax = -pi;
Bx = pi;
Ay = -pi;
By = pi;
epsilon = 0.1;
Error2=10;
N = 25;
while Error2>epsilon
N = 2*N; % Number of points; N+1 discretizations 
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
Error = 10;
U = zeros(N+2,N+2); % This is the initial guess matrix.
C = zeros(N+2,N+2); % This is the matrix that will copy the previous time frame matrix 
for i = 1:(N+2)
        C(1,i) = cos(pi*(x(i)-Ax))*cosh(Bx-x(i));
        C((N+2),i) = ((x(i)-Ax))^2*sin(pi*(x(i)-Ax)/(2*(Bx-Ax)));
     

end
    while Error>epsilon
        for i = 2:(N+1)
            for j = 2:(N+1)
                F(i,j) = cos(pi/2.*(2.*((x(j)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i)-Ay)/(By-Ay)));
                C(i,j) = 0.25*(C(i-1,j)+C(i,j-1)+U(i,j+1)+U(i+1,j)-(F(i,j)*h^2));
            end
            F(i,1) = cos(pi/2.*(2.*((x(1)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i))-Ay)/(By-Ay));
            C(i,1) = 0.25*(C(i-1,1)+C(i,2)+U(i,2)+U(i+1,1)-(F(i,1)*h^2));
            F(i,N+2) = cos(pi/2.*(2.*((x(N+2)-Ax)/(Bx-Ax))+1)).*sin(pi.*((y(i))-Ay)/(By-Ay));
            C(i,(N+2)) = 0.25*(C(i-1,N+2)+C(i,N+1)+U(i,N+1)+U(i+1,N+2)-(F(i,N+2)*h^2));
            
        end
        Err = max((abs(C-U))./max(U));
        Error = max(Err);
        P = U;
        U = C;
        counter = counter+1;
    end
    Err2 = max(abs(P-U)./max(U));
    Error2 = max(Err2);
end
    figure(1)
    surf(x,y,U)
    title(['Solution of Poisson Equation Using Gauss-Seidel Method - Iteration ',num2str(counter)]);
    xlabel('X','Fontsize',14)
    ylabel('Y','Fontsize',14)
    zlabel('U','Fontsize',14)
    Err2 = max(abs(P-U)./max(U));
    Error2 = max(Err2);

% PART B - F = 0

for i = 1:(N+2)
        C1(1,i) = cos(pi*(x(i)-Ax))*cosh(Bx-x(i));
        C1((N+2),i) = ((x(i)-Ax))^2*sin(pi*(x(i)-Ax)/(2*(Bx-Ax)));
end
Error1 = 10;
counter1 = 0;
U1 = zeros(N+2,N+2);
while Error1>epsilon
    for i = 2:(N+1)
       for j = 2:(N+1)
             C1(i,j) = 0.25*(C1(i-1,j)+C1(i,j-1)+U1(i,j+1)+U1(i+1,j)); 
       end
        C1(i,1) = 0.25*(C1(i-1,1)+C1(i,2)+U1(i,2)+U1(i+1,1));
        C1(i,(N+2)) = 0.25*(C1(i-1,N+2)+C1(i,N+1)+U1(i,N+1)+U1(i+1,N+2));
    end
Err1 = max(abs(C1-U1))./max(U1);
Error1 = max(Err1);
U1 = C1;
counter1 = counter1+1;
end
figure(2)
surf(x,y,U1)
title(['Solution of Poisson Equation Using Gauss-Seidel Method (F=0) - Iteration ',num2str(counter1)]);
xlabel('X','Fontsize',14)
ylabel('Y','Fontsize',14)
zlabel('U','Fontsize',14)