clear; clc; close all;

%% ===================== BASIC GEOMETRY & PARAMETERS =====================
R = 0.15;        % Shell radius (m)
L = 0.60;        % Length (m)
h = 0.003;       % Thickness (m)
m = 1;           % Axial mode number
n = 6;           % Circumferential mode number
p = 1;           % FGM power index
Tm = 300;        % Inner temperature (K)
Tc = 600;        % Outer temperature (K)
B0 = 0.1;        % Magnetic field (Tesla)

a0_values = [0.001 0.003 0.004]; % Initial amplitudes

%% ===================== SYMBOLIC VARIABLES =====================
syms z x theta T

%% ===================== TEMPERATURE-DEPENDENT MATERIAL LAWS =====================
Ec_sym = 406.783e9 - 22.61e6*T;                 % SiC modulus
Em_sym = 65.144e9 + 73.432e6*T - 0.1618e6*T^2;  % Al modulus

alpha_c_sym = 3e-6 + 3e-9*T - 6e-13*T^2;         % SiC CTE
alpha_m_sym = 2e-5 + 6e-9*T + 3e-12*T^2 + 1e-14*T^3; % Al CTE

rho_c = 3210; rho_m = 2700;                     % densities
mu_c  = 0.2;  mu_m = 0.3;                        % Poisson ratios
phi_c = 0;   phi_m = 3.63e7;                     % electrical conductivity

kc = 183.78; km = 235;                           % thermal conductivities
kcm = kc - km;                                   % difference

%% ===================== NONLINEAR TEMPERATURE FIELD (Eq.12) =====================
beta_sym = 0; eta_sym = 0;

for i = 0:5
    beta_sym = beta_sym + (-1)^i * (kcm^i)/((i*p + 1) * km^i);
    eta_sym  = eta_sym  + (-1)^i * (kcm^i)/((i*p + 1) * km^i) ...
               * ((2*z + h)/(2*h))^(i*p+1);
end

eta_sym = eta_sym / beta_sym;

Tz = simplify((Tc - Tm)*eta_sym + Tm);  % Temperature distribution T(z)

%% ===================== SUBSTITUTE T(z) INTO MATERIAL PROPERTIES =====================
Ec_z = simplify(subs(Ec_sym, T, Tz));
Em_z = simplify(subs(Em_sym, T, Tz));

alpha_c_z = simplify(subs(alpha_c_sym, T, Tz));
alpha_m_z = simplify(subs(alpha_m_sym, T, Tz));

%% ===================== VOLUME FRACTION =====================
Vc = (0.5 + z/h)^p;

%% ===================== EFFECTIVE PROPERTIES =====================
E_sym     = simplify((Ec_z - Em_z).*Vc + Em_z);
alpha_sym = simplify((alpha_c_z - alpha_m_z).*Vc + alpha_m_z);
rho_sym   = simplify((rho_c - rho_m).*Vc + rho_m);
phi_sym   = simplify((phi_c - phi_m).*Vc + phi_m);
mu_sym    = simplify((mu_c - mu_m).*Vc + mu_m);

%% ===================== NEUTRAL SURFACE z0 =====================
z0 = double( vpaintegral(z.*E_sym, z, -h/2, h/2) / ...
             vpaintegral(E_sym, z, -h/2, h/2) );

%% ===================== STIFFNESS MATRICES =====================
Psi_sym = E_sym./(1 - mu_sym.^2);

A11 = double(vpaintegral(Psi_sym, z, -h/2, h/2));    % A11
A22 = A11;
A12 = double(vpaintegral(mu_sym.*E_sym./(1 - mu_sym.^2), z, -h/2, h/2));
A66 = double(vpaintegral(E_sym./(2*(1 + mu_sym)), z, -h/2, h/2));

D11 = double(vpaintegral((z - z0).^2 .* E_sym./(1 - mu_sym.^2), z, -h/2, h/2));
D22 = D11;
D12 = double(vpaintegral((z - z0).^2 .* mu_sym.*E_sym./(1 - mu_sym.^2), z, -h/2, h/2));
D66 = double(vpaintegral((z - z0).^2 .* E_sym./(2*(1 + mu_sym)), z, -h/2, h/2));

%% ===================== THERMAL RESULTANTS =====================
NT = double(vpaintegral(Psi_sym.*alpha_sym.*(Tz - Tm), z, -h/2, h/2));
MT = double(vpaintegral((z-z0).*Psi_sym.*alpha_sym.*(Tz - Tm), z, -h/2, h/2));

%% ===================== MASS & MAGNETIC TERMS =====================
I0 = double(vpaintegral(rho_sym, z, -h/2, h/2));
K2 = double(vpaintegral(phi_sym*(z - z0)^2, z, -h/2, h/2));

%% ===================== LINEAR NATURAL FREQUENCY =====================
mpL = m*pi/L;    % axial wavenumber

n1_tilde = (B0^2 * K2 / I0) * (mpL^2 + n^2 / R^2);

omega0_sq = (D11/I0)*mpL^4 ...
          + (2*n^2*D12/(I0*R^2))*mpL^2 ...
          + (D22*n^4)/(I0*R^4) ...
          + (4*n^2*D66/(I0*R^2))*mpL^2 ...
          + (A22/(I0*R^2)) ...
          - (NT/I0)*(mpL^2 + n^2/R^2);

omega0 = sqrt(omega0_sq);

%% ===================== GEOMETRIC INTEGRALS L1–L4 =====================
L1 = double(vpaintegral(vpaintegral(sin(mpL*x)^2*cos(n*theta)^2, theta, 0, 2*pi), x, 0, L));
L2 = double(vpaintegral(vpaintegral(sin(mpL*x)^2*cos(mpL*x)^2*cos(n*theta)^4, theta, 0, 2*pi), x, 0, L));
L3 = double(vpaintegral(vpaintegral(sin(mpL*x)^2*cos(mpL*x)^2*sin(n*theta)^2*cos(n*theta)^2, theta, 0, 2*pi), x, 0, L));
L4 = double(vpaintegral(vpaintegral(sin(mpL*x)^4*sin(n*theta)^2*cos(n*theta)^2, theta, 0, 2*pi), x, 0, L));

%% ===================== NONLINEAR COEFFICIENT n2_tilde =====================
n2_tilde = (3*A11*mpL^4*L2)/(2*I0*L1) ...
         - (2*n^2*A12*mpL^2*L3)/(I0*R^2*L1) ...
         + (n^2*A12*mpL^2*L4)/(2*I0*R^2*L1) ...
         + (3*n^4*A22*L4)/(2*I0*R^4*L1) ...
         + (n^2*A66*L2)/(I0*R^2*L1);

%% ===================== TIME RESPONSE =====================
t = linspace(0, 40, 400);
figure; hold on;

for a0 = a0_values
    C = a0^(-2) - (3*n2_tilde)/(8*omega0^2);      % integration constant

    a_t = ( C*exp(n1_tilde*t) + ...
           (3*n2_tilde)/(8*omega0^2) ).^(-1/2);  % amplitude vs time

    omega_t = omega0 ...
        - (n1_tilde^2)/(8*omega0) ...
        + (3*n2_tilde*a_t.^2)/(8*omega0) ...
        - (15*n2_tilde^2*a_t.^4)/(256*omega0^3);

    plot(t, omega_t, 'LineWidth', 2);
end

xlabel('Time t (s)');
ylabel('\omega (rad/s)');
title('Nonlinear Natural Frequency vs Time (p = 1, n = 6)');
legend('a_0 = 0.001','a_0 = 0.003','a_0 = 0.004');
grid on; axis tight;

% for a0 = a0_values
% 
%     % ====== Integration constant C ======
%     C = a0^(-2) - (3*n2_tilde)/(8*omega0^2);
% 
%     % ====== Compute amplitude at t = 0 ======
%     denom0 = C + (3*n2_tilde)/(8*omega0^2);
% 
%     if denom0 <= 0
%         warning('a(0) becomes complex for a0 = %.6f', a0);
%     end
% 
%     a0_t = denom0^(-1/2);   % a(0)
% 
%     % ====== Compute omega(0) ======
%     omega_t0 = omega0 ...
%         - (n1_tilde^2)/(8*omega0) ...
%         + (3*n2_tilde*a0_t.^2)/(8*omega0) ...
%         - (15*n2_tilde^2*a0_t.^4)/(256*omega0^3);
% 
%     fprintf('For a0 = %.6f   a(0) = %.6f   omega(0) = %.6f rad/s\n', ...
%             a0, a0_t, omega_t0);
% 
%     % ====== Full time history for plotting ======
%     a_t = ( C*exp(n1_tilde*t) + (3*n2_tilde)/(8*omega0^2) ).^(-1/2);
% 
%     omega_t = omega0 ...
%         - (n1_tilde^2)/(8*omega0) ...
%         + (3*n2_tilde*a_t.^2)/(8*omega0) ...
%         - (15*n2_tilde^2*a_t.^4)/(256*omega0^3);
% 
%     plot(t, omega_t, 'LineWidth', 2);
% 
% end


%% ===================== END OF SCRIPT =====================================
