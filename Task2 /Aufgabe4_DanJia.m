clear;

% Define the range in X direction
X =20.; % Length in X direction

% Parameters initialized in finite difference method
n = 420;  % Number of inner gridding nodes, The total gridding nodes are n+2 including boundaries
dh = X/(n+1); % The length of step sizes

% Initial step sizes in x direction
for i = 1:n+2
    x(i) = (i-1)*dh - 10;
end

% Initial boundary conditions in X direction
for i = 1:n+2
    u(i) = 0.;
end

% Create kinetic energy operator/matrix
KenOper = (1/dh^2)*full(gallery('tridiag',n,-1,2,-1));

% Create potential energy operator/matrix, choose different potential by
% switching true ot false
Harmonic_Oscillator_potential = true;   % The potential for Aufgabe (a)
%Harmonic_Oscillator_potential = false;
%finite_square_well = true;             % The potential for Aufgabe (b)
finite_square_well = false;
%finite_square_barrier = true;          % The potential for Aufgabe (c)
finite_square_barrier = false;

if Harmonic_Oscillator_potential
    for i = 2:n+1
    PotOper(i-1,i-1) = x(i)^2;
    end
end

if finite_square_well
    PotOper = eye(n);
    Temp = zeros(n);
    for i = 2:n+1
        if x(i) >= -1 & x(i) <= 1
         Temp(i-1,i-1)= 1;
        end
    end
    PotOper = PotOper - Temp;
end

if finite_square_barrier
    PotOper = zeros(n);
    for i = 2:n+1
        if x(i) >= -1 & x(i) <= 1
         PotOper(i-1,i-1)= 1;
        end
    end   
end

% Calculate Hamiltonian, Eigenvectors and Eigenvalues by using Matlab's built in eig function
Hamiltonian = KenOper + PotOper;
[Eigenvectors, Eigenvalues] = eig(Hamiltonian);

% Plot Eigenfunctions
% It finds the lowest eigenvalues (corresponding to the ground state and plots each associated eigenfunction.
X = linspace(-10,10,n);
% If want to place at eigenvalues + Eigenvalues(1,1)*ones(n,1);
plot(X, Eigenvectors(:,1), 'g', 'LineWidth', 2)
legend('Ground State')
xlabel('Position')
ylabel('Energy')
title('The ground energy level of the lowest eigenvalue')


        
            







