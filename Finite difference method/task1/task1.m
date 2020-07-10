clear;

% Define the range in X and Y directions
X = 1.; % Length in X direction
Y = 1.; % Length in Y direction

% Parameters initialized in finite difference method
n = 100;  % Number of inner gridding nodes, The total gridding nodes are n+2 including boundaries
dh = X/(n+1); % The length of step sizes, it is equal to Y/(n+1) 

% Initial boundary conditions in X direction
for i = 1:n+2
    x(i) = (i-1)*dh;
    u(i,1) = 0.;
    u(i,n+2) = 0.;
end

% Initial boundary conditions in Y direction
for j = 1:n+2
    y(j) = (j-1)*dh;
    u(1,j) = 0.;
    u(n+2,j) = 0.;
end

% Implementation of the finite difference method
for j = 2:n+1 % Y loop
    for i = 2:n+1 % X Loop
        u(i,j) = (dh^2 * (x(i)*(1 - x(i)) + y(j)*(1-y(j))) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))/4;
    end
end

% Illustration of calculated u(i,j) using finite difference method
figure(1)
[X,Y] = meshgrid(dh:dh:1-dh);
Z = u(2:n+1,2:n+1);
surf(Z)
colorbar
xlabel("X")
ylabel("Y")
zlabel("Intensity (a.u.)")
title('Numerical solution of differential equation')

% Calculate the Fehler between numerical and analytic method at different
% positions x and y, difine the anayltic solution as Ua = 0.5*x(i)*y(j)*((x(i) - 1)*(y(j) - 1)
for j = 2:n+1 % Y loop
    for i = 2:n+1 % X Loop
        Temp(i,j) = abs(u(i,j) - (0.5*x(i)*y(j)*((x(i) - 1)*(y(j) - 1))));
    end
end

% Calculate the maximum value of Fehler
Fehler = max(max(Temp()))
Convergenceorder = log(Fehler)/log(dh)

     
        








        
        
        
        
        






