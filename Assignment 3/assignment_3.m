% 2D EM Scattering - TM Excitation (Method of Moments)
% Works for general shapes (Lines, Rectangles, L-shapes, etc.)
clear; clc; close all;

fprintf('--- 2D TM Scattering Solver (PEC) ---\n');

% --- 1. Physical Constants & Inputs ---
f = input('Enter frequency (Hz) [e.g., 300e6]: ');
phi_inc_deg = input('Enter incident angle (deg) [0 = +x direction]: ');

E0 = 1.0;   
mu0 = 4*pi*1e-7;
eps0 = 8.854e-12;
eta0 = sqrt(mu0/eps0); % Intrinsic impedance (~377 ohms)
c = 1/sqrt(mu0*eps0);
lambda = c/f;
k0 = 2*pi/lambda;
omega = 2*pi*f;
phi_inc = deg2rad(phi_inc_deg);
gamma = 0.577215665; % Euler's constant

% --- 2. Geometry Setup ---
num_lines = input('Enter number of lines/edges: ');
segs_per_line = input('Enter number of segments per line: ');

all_x = []; all_y = []; all_dl = [];
total_dist = 0;
arc_mid = [];

for i = 1:num_lines
    fprintf('\nLine %d:\n', i);
    x1 = input('  Start x: '); y1 = input('  Start y: ');
    x2 = input('  End x: ');   y2 = input('  End y: ');
    
    dx = (x2 - x1) / segs_per_line;
    dy = (y2 - y1) / segs_per_line;
    dl = sqrt(dx^2 + dy^2);
    
    for j = 0:segs_per_line-1
        all_x = [all_x; x1 + (j + 0.5) * dx];
        all_y = [all_y; y1 + (j + 0.5) * dy];
        all_dl = [all_dl; dl];
        
        % Track arc length for plotting
        arc_mid = [arc_mid; total_dist + (j + 0.5) * dl];
    end
    total_dist = total_dist + (segs_per_line * dl);
end

N = length(all_x);
Z = zeros(N, N);
V = zeros(N, 1);
pref = -(omega * mu0) / 4;

% --- 3. Matrix Filling (MoM) ---
for m = 1:N
    % Incident Electric Field (TM)
    V(m) = -E0 * exp(1i * k0 * (all_x(m)*cos(phi_inc) + all_y(m)*sin(phi_inc)));
    
    for n = 1:N
        if m == n
            % Self-term with Singularity Extraction
            term = all_dl(n) * (1 - 1i * (2/pi) * (log(k0 * all_dl(n) / 4) + gamma - 1));
            Z(m,n) = pref * term;
        else
            % Off-diagonal terms (Hankel Function)
            rho_mn = sqrt((all_x(m)-all_x(n))^2 + (all_y(m)-all_y(n))^2);
            Z(m,n) = pref * all_dl(n) * besselh(0, 2, k0 * rho_mn);
        end
    end
end

% --- 4. Solution & Normalization ---
Jz = Z \ V;
H_inc_mag = E0 / eta0;
J_norm = Jz / H_inc_mag; % Normalized Current Density J/|H_inc|


%%
% --- 5. Plotting (Matching Slide Style) ---
figure('Color', 'w', 'Name', 'TM Scattering Results');

% Magnitude Plot
subplot(2,1,1);
plot(arc_mid, abs(J_norm), 'r-+', 'LineWidth', 1, 'MarkerSize', 5);
grid on;
ylabel('J_z/|H^{inc}|');
xlabel('Arc Length (m)');
title(['Current Density Magnitude (\phi^{inc} = ', num2mstr(phi_inc_deg), '^\circ)']);

% Phase Plot
subplot(2,1,2);
% Note: Using 'unwrap' can help make the phase look continuous like the slides
plot(arc_mid, unwrap(angle(J_norm)*180/pi), 'b-*', 'LineWidth', 1, 'MarkerSize', 5);
grid on;
ylabel('\angle J_z (deg)');
xlabel('Arc Length (m)');
title(['Current Density Phase (f = ', num2str(f/1e6), ' MHz)']);
%%%%%%%%%%%%%%%%%%%%%%%%%




for m = 1:N
    V(m) = -E0 * exp(1i * k0 * (all_x(m)*cos(phi_inc) + all_y(m)*sin(phi_inc)));
    for n = 1:N
        if m == n
            term = all_dl(n) * (1 - 1i * (2/pi) * (log(k0 * all_dl(n) / 4) + gamma - 1));
            Z(m,n) = pref * term;
        else
            rho_mn = sqrt((all_x(m)-all_x(n))^2 + (all_y(m)-all_y(n))^2);
            Z(m,n) = pref * all_dl(n) * besselh(0, 2, k0 * rho_mn);
        end
    end
end
Jz = Z \ V;



