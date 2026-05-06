%% Assignment #4b - Part 1: Imbalance Axial Ratio Analysis
% Quadrature Coupler connected to a dual-feed Circularly Polarized Antenna
% Analyzes the effect of magnitude and phase imbalance on Axial Ratio (AR)
%
% The Axial Ratio formula for a CP antenna with two orthogonal components:
%
%        1 + sqrt(1 - 4*(Ex*Ey*sin(Dphi)/(Ex^2+Ey^2))^2)
% AR  = --------------------------------------------------
%        1 - sqrt(1 - 4*(Ex*Ey*sin(Dphi)/(Ex^2+Ey^2))^2)
%
% Where:
%   Ex, Ey : amplitudes of the two orthogonal E-field components
%   Dphi   : phase difference between them (phase error from 90 deg)
%
% For perfect circular polarization: Ex = Ey, Dphi = 90 deg => AR = 1 (0 dB)

close all; clear; clc;

%% =========================================================
%  HELPER FUNCTION: Compute Axial Ratio
%  =========================================================
% AR in linear, then converted to dB
AR_func = @(Ex, Ey, Dphi_deg) compute_AR_dB(Ex, Ey, Dphi_deg);

%% =========================================================
%  PLOT 1: AR vs Magnitude Imbalance (Phase fixed at 90 deg = perfect)
%  =========================================================
% Reference amplitude Ey = 1, vary Ex from 1 to ~2.5 (0 to 8 dB imbalance)
mag_imbalance_dB = linspace(0, 8, 500);   % dB
Ey = 1;
Ex_vals = Ey * 10.^(mag_imbalance_dB/20); % convert dB to linear ratio
Dphi_perfect = 90;                          % perfect quadrature phase

AR_mag = arrayfun(@(Ex) AR_func(Ex, Ey, Dphi_perfect), Ex_vals);

figure('Name','AR vs Magnitude Imbalance','Position',[100 100 700 500]);
plot(mag_imbalance_dB, AR_mag, 'b-', 'LineWidth', 2);
xlabel('Magnitude Imbalance (dB)', 'FontSize', 13);
ylabel('Axial Ratio (dB)', 'FontSize', 13);
title('AR vs Magnitude Imbalance at 90° Phase Difference', 'FontSize', 14);
grid on; grid minor;
xlim([0 8]); ylim([0 9]);
set(gca,'FontSize',12);
saveas(gcf,'AR_vs_Magnitude.png');

%% =========================================================
%  PLOT 2: AR vs Phase Imbalance (Magnitude balanced: Ex = Ey)
%  =========================================================
phase_deg = linspace(0, 180, 500);   % phase difference in degrees
Ex_bal = 1; Ey_bal = 1;             % equal amplitudes

AR_phase = arrayfun(@(ph) AR_func(Ex_bal, Ey_bal, ph), phase_deg);

figure('Name','AR vs Phase Imbalance','Position',[100 100 700 500]);
plot(phase_deg, AR_phase, 'r-', 'LineWidth', 2);
xlabel('Phase Imbalance (degrees)', 'FontSize', 13);
ylabel('Axial Ratio (dB)', 'FontSize', 13);
title('AR vs Phase Imbalance at 0 dB Magnitude Difference', 'FontSize', 14);
grid on; grid minor;
xlim([0 180]); ylim([0 45]);
set(gca,'FontSize',12);
saveas(gcf,'AR_vs_Phase.png');

%% =========================================================
%  PLOT 3: 3D Surface - AR vs Both Magnitude and Phase Imbalance
%  =========================================================
mag_imb_dB_3d = linspace(0, 8, 80);
phase_imb_3d  = linspace(0, 180, 80);

[MAG3D, PH3D] = meshgrid(mag_imb_dB_3d, phase_imb_3d);
AR_3D = zeros(size(MAG3D));

for i = 1:size(MAG3D,1)
    for j = 1:size(MAG3D,2)
        Ex_ij = 10^(MAG3D(i,j)/20);
        Ey_ij = 1;
        AR_3D(i,j) = AR_func(Ex_ij, Ey_ij, PH3D(i,j));
    end
end

figure('Name','3D AR Surface','Position',[100 100 900 600]);
surf(MAG3D, PH3D, AR_3D, 'EdgeColor','none');
colormap(jet); colorbar;
xlabel('Magnitude Imbalance (dB)', 'FontSize', 12);
ylabel('Phase Imbalance (degrees)', 'FontSize', 12);
zlabel('Axial Ratio (dB)', 'FontSize', 12);
title('3D Surface: AR vs Magnitude and Phase Imbalance', 'FontSize', 14);
view(45, 30);
set(gca,'FontSize',11);
saveas(gcf,'AR_3D_Surface.png');

%% =========================================================
%  PLOT 4: AR Contour Map (top-down view) - useful complement
%  =========================================================
figure('Name','AR Contour Map','Position',[100 100 750 550]);
contourf(MAG3D, PH3D, AR_3D, 20);
colormap(jet); colorbar;
xlabel('Magnitude Imbalance (dB)', 'FontSize', 13);
ylabel('Phase Imbalance (degrees)', 'FontSize', 13);
title('Contour Map: AR (dB) vs Magnitude & Phase Imbalance', 'FontSize', 14);
hold on;
contour(MAG3D, PH3D, AR_3D, [3 3], 'k--', 'LineWidth', 2);  % 3dB AR line
legend('AR contours','3 dB AR limit','Location','best');
set(gca,'FontSize',12);
saveas(gcf,'AR_Contour.png');

%% =========================================================
%  PLOT 5: AR vs Phase for several magnitude imbalances (family of curves)
%  =========================================================
mag_levels_dB = [0, 1, 2, 3, 4, 6];   % dB
phase_sweep = linspace(0, 180, 500);
colors_p = lines(length(mag_levels_dB));

figure('Name','AR vs Phase - Family of Curves','Position',[100 100 750 500]);
hold on;
for k = 1:length(mag_levels_dB)
    Ex_k = 10^(mag_levels_dB(k)/20);
    AR_k = arrayfun(@(ph) AR_func(Ex_k, 1, ph), phase_sweep);
    plot(phase_sweep, AR_k, 'Color', colors_p(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Mag Imb = %d dB', mag_levels_dB(k)));
end
hold off;
xlabel('Phase Imbalance (degrees)', 'FontSize', 13);
ylabel('Axial Ratio (dB)', 'FontSize', 13);
title('AR vs Phase Imbalance for Different Magnitude Imbalances', 'FontSize', 14);
legend('Location','north','FontSize',10);
grid on; grid minor;
xlim([0 180]); ylim([0 45]);
yline(3,'k--','LineWidth',1.5,'Label','3 dB limit');
set(gca,'FontSize',12);
saveas(gcf,'AR_vs_Phase_Family.png');

%% =========================================================
%  PRINT SUMMARY TABLE
%  =========================================================
fprintf('\n=== Axial Ratio Summary Table ===\n');
fprintf('%-20s %-20s %-15s\n','Mag Imb (dB)','Phase Imb (deg)','AR (dB)');
fprintf('%s\n', repmat('-',1,55));
test_mags = [0 0 0 1 2 3 4];
test_phases = [90 80 70 90 90 90 90];
for k = 1:length(test_mags)
    Ex_t = 10^(test_mags(k)/20);
    ar_t = AR_func(Ex_t, 1, test_phases(k));
    fprintf('%-20.1f %-20.1f %-15.2f\n', test_mags(k), test_phases(k), ar_t);
end

fprintf('\nAll Part 1 plots saved.\n');

%% =========================================================
%  LOCAL FUNCTION: Compute AR in dB
%  =========================================================
function AR_dB = compute_AR_dB(Ex, Ey, Dphi_deg)
    Dphi = deg2rad(Dphi_deg);
    denom = Ex^2 + Ey^2;
    inner = 1 - 4*(Ex*Ey*sin(Dphi)/denom)^2;
    inner = max(inner, 0);   % clamp to avoid numerical issues near perfect CP
    sq = sqrt(inner);
    num   = 1 + sq;
    den   = 1 - sq;
    if den <= 0
        AR_linear = Inf;
    else
        AR_linear = num/den;
    end
    % Convert to dB (AR >= 1 always; AR=1 => 0 dB = perfect CP)
    AR_dB = 20*log10(AR_linear);
    % Cap at 45 dB for display
    AR_dB = min(AR_dB, 45);
end
