% 2D EM Scattering - TM Excitation (PEC)
% Full Solver: Surface Current + Near Field Distribution
clear; clc; close all;

% --- 1. Physical Constants & Inputs ---
f = input('Enter frequency (Hz) [e.g., 300e6]: ');
phi_inc_deg = input('Enter incident angle (deg) [0 = +x direction]: ');

E0 = 1.0; 
mu0 = 4*pi*1e-7; eps0 = 8.854e-12;
eta0 = sqrt(mu0/eps0); c = 1/sqrt(mu0*eps0);
lambda = c/f; k0 = 2*pi/lambda; omega = 2*pi*f;
phi_inc = deg2rad(phi_inc_deg);
gamma = 0.577215665; 

% --- 2. Geometry Setup ---
num_lines = input('Enter number of lines/edges: ');
segs_per_line = input('Enter number of segments per line: ');

all_x = []; all_y = []; all_dl = [];
for i = 1:num_lines
    fprintf('Line %d:\n', i);
    x1 = input('  Start x: '); y1 = input('  Start y: ');
    x2 = input('  End x: ');   y2 = input('  End y: ');
    dx = (x2 - x1) / segs_per_line;
    dy = (y2 - y1) / segs_per_line;
    dl = sqrt(dx^2 + dy^2);
    for j = 0:segs_per_line-1
        all_x = [all_x; x1 + (j + 0.5) * dx];
        all_y = [all_y; y1 + (j + 0.5) * dy];
        all_dl = [all_dl; dl];
    end
end

% --- 3. MoM Matrix Solution (Surface Currents) ---
N = length(all_x);
Z = zeros(N, N);
V = -E0 * exp(1i * k0 * (all_x*cos(phi_inc) + all_y*sin(phi_inc)));
pref = -(omega * mu0) / 4;

for m = 1:N
    for n = 1:N
        if m == n
            term = all_dl(n) * (1 - 1i * (2/pi) * (log(k0 * all_dl(n) / 4) + gamma - 1));
            Z(m,n) = pref * term;
        else
            dist = sqrt((all_x(m)-all_x(n))^2 + (all_y(m)-all_y(n))^2);
            Z(m,n) = pref * all_dl(n) * besselh(0, 2, k0 * dist);
        end
    end
end
Jz = Z \ V;

% --- 4. FIXED: Optimized Field Calculation ---
fprintf('\nComputing field grid (Vectorized)...\n');
res = 150; % Grid resolution (e.g., 150x150)
margin = 1.5 * lambda; 
x_vec = linspace(min(all_x)-margin, max(all_x)+margin, res);
y_vec = linspace(min(all_y)-margin, max(all_y)+margin, res);
[X, Y] = meshgrid(x_vec, y_vec);

% Calculate Incident Field for the whole grid at once
Ei = E0 * exp(1i * k0 * (X*cos(phi_inc) + Y*sin(phi_inc)));

% Calculate Scattered Field (Es)
Es = zeros(size(X));
for n = 1:N
    % Distance from current segment n to every grid point (X,Y)
    R = sqrt((X - all_x(n)).^2 + (Y - all_y(n)).^2);
    
    % Mask points too close to the segment to avoid infinity
    % Set to a small value so besselh returns a finite (though large) number
    R(R < all_dl(n)/10) = all_dl(n)/10; 
    
    % Add contribution of segment n
    Es = Es + pref * Jz(n) * all_dl(n) * besselh(0, 2, k0 * R);
end

Etot = Ei + Es;
Etot_norm = Etot / E0;

% --- 5. Final Plotting (Matching Slide Visuals) ---
figure('Color', 'w', 'Position', [100, 100, 1000, 450]);

% Plot Real Part: Shows the Wave Pattern (Crests and Troughs)

% Plot Magnitude: Shows the Field Intensity/Shadowing
pcolor(X, Y, abs(Etot_norm));
shading interp; colormap('jet'); colorbar;
hold on;
plot(all_x, all_y, 'k-', 'LineWidth', 2); % Draw PEC object
title('|E_z^{total}| (Field Intensity)');
xlabel('x (m)'); ylabel('y (m)');
axis equal tight;