%% Assignment #4b - Master Script
% Passives: Quadrature Coupler & Phased Arrays - Beam Squint
%
% Run this script to execute both parts of the assignment.
% All figures will be saved as PNG files in the current directory.
%
% REQUIREMENTS: MATLAB (no additional toolboxes required)
%
% -------------------------------------------------------
% PART 1: Imbalance Axial Ratio Analysis
%   - Quadrature coupler feeding a dual-feed CP antenna
%   - Effect of magnitude imbalance on AR (Phase = 90°)
%   - Effect of phase imbalance on AR (Magnitude = balanced)
%   - 3D surface of AR vs both imbalances
%   - Contour map for 3 dB AR degradation boundary
%   - Family of AR curves for multiple magnitude imbalance levels
%
% PART 2: Beam Squint Analysis
%   - Radiation patterns: Phase Shifters vs TTD at 3 frequencies
%   - Beam squint vs Array size N
%   - Beam squint vs Bandwidth (multiple center frequencies)
%   - Beam squint vs Bandwidth (multiple scan angles)
%   - Beam squint vs Frequency (multiple fc values)
% -------------------------------------------------------

fprintf('========================================\n');
fprintf(' Assignment 4b - Running Part 1...\n');
fprintf('========================================\n');
run('Assignment4b_Part1_AxialRatio.m');

fprintf('\n========================================\n');
fprintf(' Assignment 4b - Running Part 2...\n');
fprintf('========================================\n');
run('Assignment4b_Part2_BeamSquint.m');

fprintf('\n========================================\n');
fprintf(' All plots generated successfully!\n');
fprintf('========================================\n');

% List generated files
fprintf('\nGenerated output files:\n');
files = dir('*.png');
for k = 1:length(files)
    fprintf('  - %s\n', files(k).name);
end
