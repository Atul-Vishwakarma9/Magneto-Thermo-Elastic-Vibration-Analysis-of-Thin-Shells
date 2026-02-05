clc;
clear;
close all;

% Geometry
L = 1.0;
R = 0.5;
h = 0.02;

% Mode numbers
m = 1;
n = 6;

% Magnetic and thermal parameters
B0 = 0.2;          % Tesla
DeltaTc = 300;    % Ceramic temperature
DeltaTm = 50;     % Metal temperature

% Initial amplitude
a0 = 0.003;

% Material law
gradation = "power";  % power / sigmoid / hyperbolic

% Compute material properties
props = material_properties(h, DeltaTc, DeltaTm, gradation);

% Compute stiffness matrices
[A, B, D] = stiffness_matrices(props, h);

% Galerkin reduction
params = galerkin_reduction(A, D, m, n, L, R, B0);

% Solve nonlinear ODE using multi-scale method
[t, omega] = multiscale_solution(params, a0);

% Plot results
plot_results(t, omega);
