function Phi = MAB522_solver_backward_euler(u, D, phi0, x, dt, steps, bcvec, every, containmentMethod)
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
% Output:
%   Phi:    solution matrix where Phi(n,i) = phi_i^n

%% Set up problem
N = length(phi0);                 % count number of nodes requiring calculaton
Phi = zeros(steps/every+1, N);    % initialise solution matrix
Phi(1,:) = phi0;                  % copy over initial condition
phi = phi0(:);                    % initialise current solution
%% Set up distance/width matrices from meshing su, D, phi0, x, dt, steps, bcvec, everycheme
dxw = zeros(N-1,1);
dxe = zeros(N-1,1);
delxP = zeros(N,1);

% Calculate nodular and CV widths from provided meshing scheme
% Left End
dxe(1) = x(2) - x(1);    % Left end next node (east)
delxP(1) = dxe(1)/2;     % Left end CV width

%Right end
dxw(N-1) = x(N) - x(N-1);   % Right end next node (west)
delxP(N) = dxw(N-1)/2;      % Right end CV width

%Internal deltas
for i = 2:N-1
    dxw(i-1) = x(i) - x(i-1);
    dxe(i) = x(i+1) - x(i);
    delxP(i) = (dxw(i-1)/2) + (dxe(i)/2);
end

%% Construct coeficcient matrix (A matrix)
diags = zeros(N,3); % Empty 3 wide matrix (to be turned into tri-diagonal system)
for j = 1:N
   if j == 1            % Start of linear surface
        diags(j,1) = -( (dt/delxP(j+1) ) * ( (u/2) + (D/dxw(j) ) ) );
        diags(j,2:3) = 1;
   elseif j == N        % End of linear surface
        diags(j,1:2) = 1;
        diags(j,3) = -((dt/delxP(j))*((-u/2) + (D/dxe(j-1))));
   else                 % Internal Nodes
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

%% Create modified boundary conditions as per containment strategy:
sourceConcentration = bcvec(3);
%debug
bdyChecker = zeros(steps,2)
concentrationPlot = zeros(steps,2)

oneAt = 0;
found = false;
%%Find the index to where x=1 approximately is:
for n = 1:N
    if (x(n) > 1) && ~found
        oneAt = n;
        found = true;
    end
end

if (oneAt == 0)
    %Error
end
%% Step in time
for n = 1:steps
    
    % Right hand side
    b = phi;
    time = (n*dt);  % store time seperately for easy reference
    
    %Swich methods based on containment strategy
    switch containmentMethod
        case 'sudden'
            if ( time > 1.25) & ~(sourceConcentration == 0) %If the time is past 0.25 days, and the value hasnt changed
                sourceConcentration = 0;
            end
        case 'gradual'
            if (time > 0.25)
                sourceConcentration = exp(-time + 0.25) * bcvec(3);
            end
        otherwise
            %Error
    end
    
    bdyChecker(n,:) = [time, sourceConcentration];
    
                
    % Adjust for boundary conditions
    b(1) = b(1) + ((D*sourceConcentration*dt)/(delxP(1)*bcvec(2)));
    b(N) = b(N) + ((D*bcvec(6)*dt)/(delxP(N)*bcvec(5)));
    
    % Compute solution at next time
    phi = U \ (L \ b(p));
    concentrationPlot(n,:) = [time, phi(oneAt)];
    
    
    % Copy solution to output if required
    if mod(n, every) == 0
        Phi(n/every + 1, :) = phi;
    end
    
end

bdyChecker
phi;

%% Plot concentration
set(figure, 'Position', get(0,'Screensize')); 
fontsize = 20;
markersize = 15;
linewidth = 1.5;

% Initial condition first (to get the legend correct)
hold on
plot(concentrationPlot(:,:), 'k.', 'MarkerSize', markersize);

% Labels, legend and title
set(xlabel('x'), 'FontSize', fontsize);
set(ylabel('y'), 'FontSize', fontsize);
set(legend('Exact', 'Location', 'Best'), 'FontSize', fontsize);
title(['Exercise 2.5 (', upper(char(method)), '). ', errorstring], 'FontSize', fontsize, 'Interpreter', 'none');
set(gca, 'FontSize', fontsize);