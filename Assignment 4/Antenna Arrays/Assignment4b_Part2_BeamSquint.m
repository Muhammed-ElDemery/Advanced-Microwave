%% Assignment #4b - Part 2: Beam Squint Analysis
% Phase Shifters vs True-Time Delay (TTD) Units in a Large Wideband Array
%
% BACKGROUND:
%  A phased array uses progressive phase shifts to steer the beam to angle θ₀.
%  Phase shifters apply a FIXED phase shift (not frequency-dependent).
%  TTD units apply a time delay τ = (d/c)*sin(θ₀), which is frequency-independent.
%
%  Beam steering condition at design frequency fc:
%    ψ(fc) = k(fc)*d*sin(θ₀) = (2π*fc/c)*d*sin(θ₀)
%
%  At another frequency f, with PHASE SHIFTERS (ψ fixed):
%    sin(θ_actual) = (fc/f)*sin(θ₀)   =>  θ_actual = arcsin((fc/f)*sin(θ₀))
%    Beam Squint = θ_actual - θ₀
%
%  With TTD: beam stays at θ₀ for all frequencies (no squint)
%
%  Element spacing d = λc/2 = c/(2*fc)

close all; clear; clc;

c = 3e8;  % speed of light (m/s)

%% =========================================================
%  PART A: Radiation Patterns - Phase Shifters vs TTD
%  Compare at 3 frequencies across a wideband signal
%  =========================================================
fc = 50e9;           % center frequency (Hz) - 50 GHz (mm-wave)
BW_pct = 20;         % bandwidth %
f_low  = fc*(1 - BW_pct/200);    % 45 GHz
f_high = fc*(1 + BW_pct/200);    % 55 GHz
freqs  = [f_low, fc, f_high];    % [45, 50, 55] GHz

theta0_deg = 30;      % desired scan angle
theta0 = deg2rad(theta0_deg);

N  = 64;              % number of array elements
d  = c/(2*fc);        % half-wavelength spacing at fc

theta_scan = linspace(0, 60, 1000);   % scan range (degrees)
theta_rad  = deg2rad(theta_scan);

%-- Array factor function
AF_func = @(theta_vec, f, psi_shift) compute_AF_dB(theta_vec, f, psi_shift, N, d, c);

%-- Phase shifts set at fc for θ₀
psi_fc = 2*pi*fc/c * d * sin(theta0);   % progressive phase shift at fc

figure('Name','Phase Shifters vs TTD Patterns','Position',[50 50 1300 550]);

%--- Left: Phase Shifters
subplot(1,2,1);
hold on;
freq_labels_PS = {};
colors_f = [0 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13];
for k = 1:length(freqs)
    f_k = freqs(k);
    % Phase shifter: fixed psi_fc regardless of frequency
    AF_k = AF_func(theta_rad, f_k, psi_fc);
    
    % Find actual beam peak
    [~, idx_peak] = max(AF_k);
    theta_peak = theta_scan(idx_peak);
    squint = theta_peak - theta0_deg;
    
    plot(theta_scan, AF_k, 'Color', colors_f(k,:), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%.0f GHz (%.1f°)', f_k/1e9, theta_peak));
    freq_labels_PS{k} = sprintf('%.0f GHz (%.1f°)', f_k/1e9, theta_peak);
end
hold off;
xlabel('Angle (degrees)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('Phase Shifters: Beam Squint Effect', 'FontSize', 13, 'FontWeight','bold');
legend(freq_labels_PS, 'Location', 'northwest', 'FontSize', 9);
grid on; xlim([0 60]); ylim([-40 5]);
set(gca,'FontSize',11);
text(2, -35, 'Phase shifters are frequency-independent → beam squints', ...
    'FontSize', 9, 'Color', [0.5 0 0]);

%--- Right: TTD Units
subplot(1,2,2);
hold on;
freq_labels_TTD = {};
for k = 1:length(freqs)
    f_k = freqs(k);
    % TTD: progressive phase shift scales with frequency (time delay τ = d*sin(θ₀)/c)
    tau = d*sin(theta0)/c;
    psi_ttd = 2*pi*f_k*tau;   % frequency-dependent phase
    AF_k = AF_func(theta_rad, f_k, psi_ttd);
    
    [~, idx_peak] = max(AF_k);
    theta_peak = theta_scan(idx_peak);
    
    plot(theta_scan, AF_k, 'Color', colors_f(k,:), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%.0f GHz (%.1f°)', f_k/1e9, theta_peak));
    freq_labels_TTD{k} = sprintf('%.0f GHz (%.1f°)', f_k/1e9, theta_peak);
end
hold off;
xlabel('Angle (degrees)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('True-Time Delay Units: No Beam Squint', 'FontSize', 13, 'FontWeight','bold');
legend(freq_labels_TTD, 'Location', 'northwest', 'FontSize', 9);
grid on; xlim([0 60]); ylim([-40 5]);
set(gca,'FontSize',11);
text(2, -35, 'TTD gives equal time delay → beam stays at θ₀ = 30°', ...
    'FontSize', 9, 'Color', [0 0.5 0]);

sgtitle('Phase Shift vs TTD: Radiation Patterns at Multiple Frequencies', ...
    'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf,'PhaseShifter_vs_TTD_Patterns.png');

%% =========================================================
%  PART B: Beam Squint vs Array Size (N)
%  =========================================================
N_vals = 10:5:120;
theta0_bs = 30;
BW_bs = 20;   % % bandwidth
fc_bs = 50e9;
% Beam squint formula (analytical):
% sin(theta_squint) = (fc/f)*sin(theta0) - sin(theta0) at band edge
% => squint = arcsin((fc/f_edge)*sin(theta0)) - theta0
% Note: squint is independent of N (for uniform linear array)

f_edge = fc_bs*(1 + BW_bs/200);   % upper band edge gives negative squint
squint_analytic = rad2deg(asin((fc_bs/f_edge)*sind(theta0_bs))) - theta0_bs;

squint_N = zeros(size(N_vals));
for idx = 1:length(N_vals)
    squint_N(idx) = compute_squint_numerical(N_vals(idx), fc_bs, f_edge, theta0_bs, c);
end

figure('Name','Beam Squint vs Array Size','Position',[100 100 700 500]);
plot(N_vals, squint_N, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'DisplayName', 'Beam Squint Angle');
yline(squint_analytic, 'r--', 'LineWidth', 1.5, 'Label', 'Analytical');
xlabel('Number of Array Elements (N)', 'FontSize', 13);
ylabel('Beam Squint Angle (degrees)', 'FontSize', 13);
title(sprintf('Beam Squint vs Array Size\n(f_c = %.0f GHz, BW = %d%%, θ_0 = %d°)', ...
    fc_bs/1e9, BW_bs, theta0_bs), 'FontSize', 13);
legend('Location','best','FontSize',11);
grid on; grid minor;
xlim([10 120]);
set(gca,'FontSize',12);
fprintf('\nNote: Beam squint is nearly independent of array size N.\n');
fprintf('Analytical squint at upper BW edge: %.2f degrees\n', squint_analytic);
saveas(gcf,'BeamSquint_vs_N.png');

%% =========================================================
%  PART C: Beam Squint vs Bandwidth
%  =========================================================
BW_range = linspace(0, 40, 200);   % % of center frequency
N_bw = 64;
fc_bw = 50e9;
theta0_vals_bw = [30];   % degrees
freq_centers = [25e9, 50e9, 75e9];   % multiple center frequencies

figure('Name','Beam Squint vs Bandwidth','Position',[100 100 750 520]);
hold on;
colors_bw = lines(length(freq_centers));
for m = 1:length(freq_centers)
    fc_m = freq_centers(m);
    squint_bw = zeros(size(BW_range));
    for idx = 1:length(BW_range)
        bw_m = BW_range(idx);
        if bw_m == 0
            squint_bw(idx) = 0;
        else
            f_edge_m = fc_m*(1 + bw_m/200);
            arg = (fc_m/f_edge_m)*sind(30);
            if abs(arg) > 1; squint_bw(idx) = NaN; else
                squint_bw(idx) = rad2deg(asin(arg)) - 30;
            end
        end
    end
    plot(BW_range, abs(squint_bw), 'Color', colors_bw(m,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('f_c = %.0f GHz', fc_m/1e9));
end
hold off;
xlabel('Bandwidth (% of frequency)', 'FontSize', 13);
ylabel('Beam Squint Angle (degrees)', 'FontSize', 13);
title(sprintf('Beam Squint vs Bandwidth for Different Center Frequencies\n(N = %d, θ_0 = 30°)', N_bw), ...
    'FontSize', 13);
legend('Location','northwest','FontSize',11);
grid on; grid minor;
xlim([0 40]); ylim([0 9]);
set(gca,'FontSize',12);
saveas(gcf,'BeamSquint_vs_Bandwidth.png');

%% =========================================================
%  PART D: Beam Squint vs Bandwidth for Different Scan Angles
%  =========================================================
theta0_set = [-45, -30, -15, 15, 30, 45];   % degrees
BW_range2 = linspace(0, 40, 200);
fc_d = 50e9;
N_d = 64;

colors_th = lines(length(theta0_set));
linestyles = {'-','--','-.',':','-','--'};

figure('Name','Beam Squint vs BW - Scan Angles','Position',[100 100 800 550]);
hold on;
leg_str = {};
for m = 1:length(theta0_set)
    th0 = theta0_set(m);
    squint_th = zeros(size(BW_range2));
    for idx = 1:length(BW_range2)
        bw_v = BW_range2(idx);
        if bw_v == 0
            squint_th(idx) = 0;
        else
            f_edge_d = fc_d*(1 + bw_v/200);   % upper edge (squint toward endfire)
            arg = (fc_d/f_edge_d)*sind(th0);
            if abs(arg) > 1; squint_th(idx) = NaN; else
                squint_th(idx) = rad2deg(asin(arg)) - th0;
            end
        end
    end
    plot(BW_range2, squint_th, 'Color', colors_th(m,:), ...
        'LineStyle', linestyles{m}, 'LineWidth', 2, ...
        'DisplayName', sprintf('θ_0 = %d°', th0));
    leg_str{end+1} = sprintf('θ_0 = %d°', th0);
end
hold off;
xlabel('Bandwidth (% of frequency)', 'FontSize', 13);
ylabel('Beam Squint Angle (degrees)', 'FontSize', 13);
title(sprintf('Beam Squint vs Bandwidth at Different Scan Angles\n(Frequency = %.0f GHz, N = %d)', ...
    fc_d/1e9, N_d), 'FontSize', 13);
legend(leg_str, 'Location','southwest','FontSize',10,'NumColumns',2);
grid on; grid minor;
xlim([0 40]);
set(gca,'FontSize',12);
saveas(gcf,'BeamSquint_vs_BW_ScanAngles.png');

%% =========================================================
%  PART E: Beam Squint vs Frequency (fc as parameter)
%  =========================================================
freq_sweep = linspace(5e9, 100e9, 500);   % Hz
N_e = 64;
fc_vals = [25e9, 50e9, 75e9];   % design center frequencies
theta0_e_set = [-45,-30,-15,15,30,45];

for fc_idx = 1:length(fc_vals)
    fc_e = fc_vals(fc_idx);
    
    figure('Name',sprintf('Beam Squint vs Frequency fc=%.0fGHz',fc_e/1e9),...
        'Position',[100 100 800 550]);
    hold on;
    colors_e = lines(length(theta0_e_set));
    lsyles = {'-','--','-.',':','-','--'};
    
    for m = 1:length(theta0_e_set)
        th0_e = theta0_e_set(m);
        squint_freq = zeros(size(freq_sweep));
        for idx = 1:length(freq_sweep)
            f_v = freq_sweep(idx);
            arg = (fc_e/f_v)*sind(th0_e);
            arg = max(-1, min(1, arg));
            squint_freq(idx) = rad2deg(asin(arg)) - th0_e;
        end
        plot(freq_sweep/1e9, squint_freq, 'Color', colors_e(m,:), ...
            'LineStyle', lsyles{m}, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('θ_0 = %d°', th0_e));
    end
    hold off;
    xlabel('Frequency (GHz)', 'FontSize', 13);
    ylabel('Beam Squint Angle (degrees)', 'FontSize', 13);
    title(sprintf('Beam Squint vs Frequency  (F = %.0f GHz, N = %d)', ...
        fc_e/1e9, N_e), 'FontSize', 14, 'FontWeight','bold');
    legend('Location','east','FontSize',10,'NumColumns',2);
    grid on; grid minor;
    xlim([5 100]); ylim([-80 80]);
    yline(0,'k-','LineWidth',0.8);
    xline(fc_e/1e9,'k--','LineWidth',1.2,'Label',sprintf('f_c=%.0fGHz',fc_e/1e9));
    set(gca,'FontSize',12);
    saveas(gcf, sprintf('BeamSquint_vs_Freq_fc%.0fGHz.png', fc_e/1e9));
end

fprintf('\nAll Part 2 plots saved.\n');

%% =========================================================
%  LOCAL FUNCTIONS
%  =========================================================
function AF_dB = compute_AF_dB(theta_vec, f, psi_shift, N, d, c)
    k = 2*pi*f/c;
    psi = k*d*sin(theta_vec) - psi_shift;   % array phase argument
    % Array factor (uniform linear array)
    % AF = sin(N*psi/2) / (N*sin(psi/2))
    AF = sin(N*psi/2) ./ (N*sin(psi/2 + 1e-30));
    AF_dB = 20*log10(abs(AF) + 1e-10);
end

function squint = compute_squint_numerical(N, fc, f_test, theta0_deg, c)
    d = c/(2*fc);
    theta0 = deg2rad(theta0_deg);
    psi_fc = 2*pi*fc/c * d * sin(theta0);
    
    theta_scan = linspace(-90, 90, 5000);
    theta_rad = deg2rad(theta_scan);
    
    k = 2*pi*f_test/c;
    psi = k*d*sin(theta_rad) - psi_fc;
    AF = abs(sin(N*psi/2) ./ (N*sin(psi/2 + 1e-30)));
    
    % Search near theta0 for peak
    [~, idx] = max(AF);
    squint = theta_scan(idx) - theta0_deg;
end
