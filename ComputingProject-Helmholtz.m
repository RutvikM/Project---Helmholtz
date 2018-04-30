clear all,clc
Ax = -pi;
Bx = pi;
Ay = -pi;
By = pi;
N = 15; % N+1 discretization
h = (Bx-Ax)/(N+1); % Since the discretization for x and y are equal, h is used.
x = zeros(1,(N+1));
y = zeros(1,(N+1));
x(1,1) = Ax;
y(1,1) = Ay;
for a = 2:(N+1)
    x(1,a) = Ax+(a*h); % Grid being made in the for loop.
    y(1,a) = Ay+(a*h);
end
%Building the penta-diagonal Matrix
A = zeros(N,N);
a = 1;
k = 3/2;
b = -4+(k*h^2);
c = -k^2*h^2;
for l = 1:N
  A(l,l) = b;
end
for l = 1:(N-3) 
  A(l+3,l) = a;
  A(l,l+3) = a;
end
for l = 1:(N-1)
  A(l,l+1) = a;
  A(l+1,l) = a;
end
z = N-1;
d = length(z);
z = floor(z/3);
for l = 1:z
  A(1+(3*l),(3*l))= 0;
  A((3*l),1+(3*l))= 0;
end