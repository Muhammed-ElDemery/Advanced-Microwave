%% Assignment 2: Comprehensive MoM for 2D PEC Flat Strip

clear; clc; close all;
start_time = cputime;
% 1. Parameters
w = 0.5;                % Half-width (Total width = 1m)
V0 = 1;                 % Potential on strip (1 Volt)
epsilon0 = 8.854e-12;   % Permittivity of free space

% Define N values for convergence study
N_list = [10, 20, 50, 100]; 
total_charge_normalized = zeros(size(N_list));

figure(1); hold on;

for k = 1:length(N_list)
    N = N_list(k);
    dx = (2*w) / N;
    x_nodes = linspace(-w, w, N+1);
    x_mid = (x_nodes(1:end-1) + x_nodes(2:end)) / 2; % Matching points
    
    % 2. Build MoM Matrix S (Point Matching)
    % Using analytical integral of ln|xm - x'| to handle singularity
    S = zeros(N, N);
    for m = 1:N
        xm = x_mid(m);
        for n = 1:N
            a = x_nodes(n);
            b = x_nodes(n+1);
            % Integral of ln|xm - x'| from a to b
            term1 = (b-xm) * log(abs(b-xm) + eps);
            term2 = (a-xm) * log(abs(a-xm) + eps);
            integral_val = term1 - term2 - (b-a);
            S(m,n) = (-1 / (2*pi*epsilon0)) * integral_val;
        end
    end
    
    % 3. Solve for Charge Coefficients q
    q = S \ (V0 * ones(N, 1));
    
    % 4. Store Data for Plotting
    total_charge_normalized(k) = sum(q * dx) / epsilon0;
    plot(x_mid, q/epsilon0, 'DisplayName', ['N = ' num2str(N)], 'LineWidth',2);
end



% --- Task 1: Charge Distribution ---
title('Normalized Charge Distribution q(x)/\epsilon_0');
xlabel('Position x (m)'); ylabel('q(x) / \epsilon_0');
legend; grid on;
%%
% --- Task 2: Convergence Study ---
figure(2);
plot(1./N_list, total_charge_normalized, 'rs-', 'LineWidth', 2);
xlabel('1/N'); ylabel('Q_L / \epsilon_0');
title('Convergence of Total Charge (Q_L / \epsilon_0 vs 1/N)');
grid on;


elapsed_time = cputime-start_time;
fprintf('Execution time: %.6f seconds\n', elapsed_time);
fprintf("Total charge:");
total_charge_normalized(end)
% --- Task 3: Potential Error on the Strip ---
% Using the results from the highest N
x_fine = linspace(-w, w, 500);
V_on_strip = zeros(size(x_fine));
for i = 1:length(x_fine)
    xi = x_fine(i);
    for n = 1:N % Using N=100
        a = x_nodes(n); b = x_nodes(n+1);
        term1 = (b-xi) * log(abs(b-xi) + eps);
        term2 = (a-xi) * log(abs(a-xi) + eps);
        V_on_strip(i) = V_on_strip(i) + (q(n) * (-1/(2*pi*epsilon0)) * (term1 - term2 - (b-a)));
    end
end
figure(3);
plot(x_fine, V_on_strip - V0, 'b-', 'LineWidth',2); hold on;
plot(x_mid, zeros(size(x_mid)), 'ro', 'MarkerSize', 3, 'LineWidth',2);
title(['Potential Error on Strip (N = ' num2str(N) ')']);
xlabel('x (m)'); ylabel('\Phi(x) - V_0 (V)');
legend('Error', 'Matching Points'); grid on;

% --- Task 4 & 5: Potential and E-Field in Surrounding Region ---
[X, Y] = meshgrid(linspace(-2, 2, 50), linspace(-2, 2, 50));
Phi = zeros(size(X)); Ex = zeros(size(X)); Ey = zeros(size(X));

for i = 1:numel(X)
    xi = X(i); yi = Y(i);
    for n = 1:N
        % Distance from midpoint of segment n to grid point (xi, yi)
        r = sqrt((xi - x_mid(n))^2 + yi^2);
        
        % Potential contribution (midpoint approximation for external points)
        pot_c = (q(n) * dx * (-1/(2*pi*epsilon0)) * log(r + eps));
        Phi(i) = Phi(i) + pot_c;
        
        % E-Field contribution E = -grad(Phi)
        E_const = (q(n) * dx) / (2*pi*epsilon0 * r^2 + eps);
        Ex(i) = Ex(i) + E_const * (xi - x_mid(n));
        Ey(i) = Ey(i) + E_const * (yi);
    end
end

figure(4);

% Filled contour for full potential
contourf(X, Y, Phi, 40);
colorbar;
hold on;

% Draw the PEC strip
%plot([-w w], [0 0], 'w', 'LineWidth', 4);

% --- Zero potential contour ---
[C,h] = contour(X, Y, Phi, [0 0], 'k', 'LineWidth', 2);

% Optional: make it thicker or different style
set(h, 'LineWidth', 2.5);

title('Potential Distribution \Phi(x,y) around the Strip');
xlabel('x (m)');
ylabel('y (m)');

% Plot Electric Field (Task 5)
figure(5);
quiver(X, Y, Ex, Ey, 2.5); hold on;
plot([-w w], [0 0], 'r', 'LineWidth', 3); % Draw the strip in red
title('Electric Field Distribution E(x,y)');
xlabel('x (m)'); ylabel('y (m)'); axis equal;

