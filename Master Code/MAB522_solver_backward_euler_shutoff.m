function Phi = MAB522_solver_backward_euler(u, D, phi0, x, dt, steps, bcvec, every)
%% 2nd order spatial accurate solution of advection-diffusion equation using backward Euler time stepping
% Solves dphi/dt + u dphi/dx = D d^2phi/dx^2
% subject to general boundary conditions
% Input:
%   u:      velocity
%   D:      diffusion coefficient
%   phi0:   initial condition vector
%   dx:     subinterval width
%   dt:     time stepsize
%   steps:  number of time steps
%   bcvec:  vector of boundary condition coefficients
%   every:  return only every ... solution
%
%   Output:
%   Phi:    solution matrix where Phi(n,i) = phi_i^n
%   shutoff shutoff method ("instant", "gradual", "none")
%

%% Set up problem
N = length(phi0);
Phi = zeros(steps/every+1, N);    % initialise solution matrix
Phi(1,:) = phi0;                  % copy over initial condition
phi = phi0(:);                    % initialise current solution
%% Set up matrix
dxw = zeros(N-1,1);
dxe = zeros(N-1,1);
delxP = zeros(N,1);

for i = 1:N
    if i == 1
        dxe(i) = x(i+1) - x(i);
        delxP(i) = dxe(i)/2;
    elseif i == N
        dxw(i-1) = x(i) - x(i-1);
        delxP(i) = dxw(i-1)/2;
    else
        dxw(i-1) = x(i) - x(i-1);
        dxe(i) = x(i+1) - x(i);
        delxP(i) = (dxw(i-1)/2) + (dxe(i)/2);
    end
end

diags = zeros(N,3);
for j = 1:N
   if j == 1
        diags(j,1) = -((dt/delxP(j+1))*((u/2) + (D/dxw(j))));
        diags(j,2) = 1;
        diags(j,3) = 1;
   elseif j == N
        diags(j,1) = 1;
        diags(j,2) = 1;
        diags(j,3) = -((dt/delxP(j))*((-u/2) + (D/dxe(j-1))));
   else
        diags(j,1) = -((dt/delxP(j))*((u/2) + (D/dxw(j-1))));
        diags(j,2) = (1 + ((D*dt)/(delxP(j)*dxw(j-1))) + ((D*dt)/(delxP(j)*dxe(j-1))));
        diags(j,3) = -((dt/delxP(j))*((-u/2) + (D/dxe(j))));
   end
end

B = spdiags(diags,-1:1,N,N);

% Adjust for boundary conditions

B(1,1) = 1 + (((-u/2) + (D/dxe(1)) + (D*bcvec(1)/bcvec(2)))*(dt/delxP(1)));
B(1,2) = -((-u/2) + (D/dxe(1)))*(dt/delxP(1));

B(N,N) = 1 + (((u/2) + (D/dxw(N-1)) + (D*bcvec(4)/bcvec(5)))*(dt/delxP(N)));
B(N,N-1) = -(((u/2) + (D/dxw(N-1)))*(dt/delxP(N)));

[L,U,p] = lu(B, 'vector');

%% Step in time
for n = 1:steps
    
    % Right hand side
    b = phi;
    
    % Adjust for boundary conditions
    b(1) = b(1) + ((D*bcvec(3)*dt)/(delxP(1)*bcvec(2)));
    b(N) = b(N) + ((D*bcvec(6)*dt)/(delxP(N)*bcvec(5)));
    
    % Compute solution at next time
    phi = U \ (L \ b(p));
    
    % Copy solution to output if required
    if mod(n, every) == 0
        Phi(n/every + 1, :) = phi;
    end

end