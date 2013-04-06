%% Solving the advection-diffusion equation
% In this exercise we are using finite differences to solve the advection-
% diffusion equation
% $$ \frac{\partial \varphi}{\partial t} + u \frac{\partial \varphi}{\partial x} = D \frac{\partial^2 \varphi}{\partial x^2} $$
% subject to general boundary conditions.

%% Clear workspace
clear;

%% Physical parameters
L = 1.1;                  % length
u = 0.5;                % velocity
D = 0.0025;             % diffusivity
c0 = 1;                 % boundary concentration

%% Boundary condition vector
H = 1000;
bcvec = [H, 1, H*c0, 0, 1, 0];
%bcvec = [H, 1, H*c0, H, 1, 0];

%% Numerical parameters
final_time = 2;                                 % final time to step to
time_steps = 500;                               % number of time steps
subintervals = 50;                              % number of subintervals in mesh
method = @MAB522_solver_backward_euler;   % numerical method

%% Visualisation parameters
num_plots = 10;             % number of plots (excluding initial condition)
subintervals_fine = 1000;   % number of subintervals in fine mesh

%% Space and time discretisation quantities for numerical solution
theta = 1;
dx = L / subintervals;                  % subinterval width
x = linspace(0, L, subintervals + 1);   % mesh
dt = final_time / time_steps;           % time step size
every = time_steps / num_plots;         % how often to output solution

%% Space and time discretisation quantities for exact solution
xfine = linspace(0, L, subintervals_fine+1);    % fine mesh
tvec = linspace(0, final_time, num_plots+1);    % vector of output times

%% Compute numerical solution
phi0 = zeros(1, subintervals + 1); phi0(1) = c0;
numerical = method(u, D, phi0, x, dt, time_steps, bcvec, every);

%% Compute exact solution on same mesh
[X,T] = meshgrid(x, tvec);
exact = exercise2_6_exact_solution(u, D, 0, c0, X, T);

%% Compute maximum error
maxerror = max(max(abs(exact(:,1:end-5) - numerical(:,1:end-5))));
if any(isnan(numerical(end,:))), maxerror = nan; end
errorstring = sprintf('Maximum error = %.3d\n', maxerror);
disp(errorstring);

%% Compute exact solution on fine mesh (for plotting)
[X,T] = meshgrid(xfine, tvec);
exact_fine = exercise2_6_exact_solution(u, D, 0, c0, X, T);

%% Plot solutions
set(figure, 'Position', get(0,'Screensize')); 
fontsize = 20;
markersize = 15;
linewidth = 1.5;

% Initial condition first (to get the legend correct)
plot(xfine, exact_fine(1,:), 'k', 'LineWidth', linewidth);
hold on
plot(x, numerical(1,:), 'k.', 'MarkerSize', markersize);

% Other solutions
plot(xfine, exact_fine(2:end, :)', 'LineWidth', linewidth);
plot(x, numerical(2:end, :)', '.', 'MarkerSize', markersize);

% Labels, legend and title
set(xlabel('x'), 'FontSize', fontsize);
set(ylabel('y'), 'FontSize', fontsize);
set(legend('Exact', 'Numerical', 'Location', 'Best'), 'FontSize', fontsize);
title(['Exercise 2.5 (', upper(char(method)), '). ', errorstring], 'FontSize', fontsize, 'Interpreter', 'none');
set(gca, 'FontSize', fontsize);