% 2D EM Scattering - TM Excitation (Method of Moments)
% Solves for Current Density Jz on a PEC surface
clear; clc; close all;

fprintf('--- 2D TM Scattering Solver (PEC) ---\n');

% --- 1. Physical Constants & Inputs ---
f = input('Enter frequency (Hz) [e.g., 300e6]: ');
phi_inc_deg = input('Enter incident angle (deg) [0 = +x direction]: ');

E0 = 1.0;
mu0 = 4*pi*1e-7;
eps0 = 8.854e-12;
eta0 = sqrt(mu0/eps0); 
c = 1/sqrt(mu0*eps0);
lambda = c/f;
k0 = 2*pi/lambda;
omega = 2*pi*f;
phi_inc = deg2rad(phi_inc_deg);
gamma = 0.5772156649; % Euler's constant

% --- 2. Geometry Setup ---
% Example: Define a simple square or custom lines
num_lines = input('Enter number of lines/edges (e.g., 4 for a square): ');
segs_per_line = input('Enter segments per line (e.g., 20): ');

all_x = []; all_y = []; all_dl = []; arc_mid = [];
total_dist = 0;

for i = 1:num_lines
    fprintf('\nLine %d:\n', i);
    x1 = input(' Start x: '); y1 = input(' Start y: ');
    x2 = input(' End x: ');   y2 = input(' End y: ');
    
    dx = (x2 - x1) / segs_per_line;
    dy = (y2 - y1) / segs_per_line;
    dl = sqrt(dx^2 + dy^2);
    
    for j = 0:segs_per_line-1
        % Midpoint of the segment (Collocation point)
        all_x = [all_x; x1 + (j + 0.5) * dx];
        all_y = [all_y; y1 + (j + 0.5) * dy];
        all_dl = [all_dl; dl];
        
        % Track cumulative arc length for plotting
        arc_mid = [arc_mid; total_dist + (j + 0.5) * dl];
    end
    total_dist = total_dist + (segs_per_line * dl);
end

N = length(all_x);
Z = zeros(N, N);
V = zeros(N, 1);

% Prefactor for the 2D Green's function (TM case)
% G = -(j/4) * H0(2)(k*rho) -> Zmn = (omega*mu/4) * dl * H0(2)(k*rho)
pref = (omega * mu0) / 4;

% --- 3. Matrix Filling (MoM) ---
for m = 1:N
    % Incident Electric Field (Traveling towards the origin)
    % Ei = E0 * exp(-j * k0 * (x*cos(phi) + y*sin(phi)))
    V(m) = E0 * exp(-1i * k0 * (all_x(m)*cos(phi_inc) + all_y(m)*sin(phi_inc)));
    
    for n = 1:N
        if m == n
            % Self-term: Integration of the Hankel function singularity
            % Using approximation: 1 - j*(2/pi)*(ln(k*dl/4e) + gamma)
            % where 1/e results in the "- 1" term inside the log
            Z(m,n) = pref * all_dl(n) * (1 - 1i * (2/pi) * (log(k0 * all_dl(n) / 4) + gamma - 1));
        else
            % Off-diagonal terms (Hankel Function)
            rho_mn = sqrt((all_x(m)-all_x(n))^2 + (all_y(m)-all_y(n))^2);
            Z(m,n) = pref * all_dl(n) * besselh(0, 2, k0 * rho_mn);
        end
    end
end

% --- 4. Solution & Normalization ---
% Solve [Z]{Jz} = {V} 
% Note: Jz here is actually the current density Jz
Jz = Z \ V;

H_inc_mag = E0 / eta0;
J_norm = Jz / H_inc_mag; % Normalized Current Density Jz/|H_inc|

% --- 5. Plotting ---
figure('Color', 'w', 'Position', [100, 100, 800, 600]);

% Magnitude Plot
subplot(2,1,1);
plot(arc_mid, abs(J_norm), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;
ylabel('|J_z / H_{inc}|');
xlabel('Arc Length (m)');
title(['Current Density Magnitude (\phi_{inc} = ', num2str(phi_inc_deg), '^\circ)']);

% Phase Plot
subplot(2,1,2);
% 1. Calculate the raw phase (in degrees)
raw_phase_deg = angle(J_norm) * 180 / pi;

% 2. Unwrap it to make it smooth (this might shift the whole thing to 300+ degrees)
unwrapped_phase = unwrap(raw_phase_deg, 170); % 170 is the jump tolerance

% 3. FIND THE OFFSET: 
% Find the index corresponding to the 'front face' (e.g., the middle of your arc)
mid_idx = round(length(unwrapped_phase) / 2); 
phase_offset = unwrapped_phase(mid_idx);

% 4. SHIFT THE PLOT so the middle is exactly 0
final_phase = unwrapped_phase - phase_offset;

% --- Plotting ---

plot(arc_mid, final_phase, 'b-', 'LineWidth', 2);
grid on;

% 5. FIX THE Y-AXIS: Explicitly set ticks and limits
ylabel('Phase \angle J_z (degrees)');
xlabel('Arc Length (m)');
title('Normalized Current Density Phase (0° at Center)');

% This ensures 0 is a main grid line and the axis looks "normal"
ylim([min(final_phase)-30, max(final_phase)+30]); 
hline = yline(0, 'r--', 'HandleVisibility', 'off'); % Optional: red line at zero