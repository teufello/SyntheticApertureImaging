%% Init
addpath(genpath('.\code'))
addpath(genpath('.\data'))
addpath(genpath('.\Field_II_ver_3_30_windows'))
% field_init();
clear all; close all;


%% Define the parameters of the transducer
N_elements = 192;
Active_elements = 64;
f0 = 7e6;
fs = 100e6;
c = 1540;
lambda = c/f0; % wavelength
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % Spacing between virtual sources = 0.23/1000 (m)
element_height = 5/1000;
focus = [0,0,30/1000];
F_recv = 1;
F_transm = -1;
pixel_size = 0.1/1000; %  pixel size laterally/axially 

%% Parameters
% Scanning
% Define how many transmissions
no_transm = 128;
% Compute the image width
image_width = (no_transm-1)*pitch;
% Compute the position of the left most line
x = -image_width/2;
% Compute how much to display per line
d_x = pitch;

% Matrix 
num_depth = 70/1000/pixel_size;
num_lateral = image_width/pixel_size;
depth_range = linspace(0.1, 70, num_depth) / 1000; % Convert
lateral_range = linspace(-image_width/2, image_width/2, num_lateral); 

%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%
%% Load data
% point phantom
load image_HRI_beamforming_point_phantom.mat
% tissue motion
% load image_HRI_beamforming_motion_lambda16.mat
% load image_HRI_beamforming_motion_lambda64.mat
rf_data  = abs(HRI);


%% PSF Analysis: Energy in Main Lobe vs. Side Lobes (10 mm x 10 mm region)

% Point target location
depth_psf = 30.3/1000;      % m
lateral_psf = 0.1/1000;    % m, steady: small value improves the positioning of the circle center; motion: choose up to 1mm

% Find center indices
[~, zc_idx] = min(abs(depth_range - depth_psf));
[~, xc_idx] = min(abs(lateral_range - lateral_psf));

% Convert 10 mm to pixels
psf_window_size = 10 / 1000; % in m
psf_win_pix = round(psf_window_size / pixel_size);
half_win = floor(psf_win_pix / 2);

% Define region around point
z_psf_range = (zc_idx - half_win):(zc_idx + half_win - 1);
x_psf_range = (xc_idx - half_win):(xc_idx + half_win - 1);

% Extract amplitude patch
psf_patch = rf_data(z_psf_range, x_psf_range);

% Create spatial grid for mask (in meters)
[xx, zz] = meshgrid(1:psf_win_pix, 1:psf_win_pix);
center = (psf_win_pix + 1) / 2;
dx = (xx - center) * pixel_size;  % lateral offset from center
dz = (zz - center) * pixel_size;  % depth offset from center
rr = sqrt(dx.^2 + dz.^2);

% Define circular mask with radius 2.5 * lambda
r_circle = 2.5 * lambda;  % in meters
mask_circle = rr <= r_circle;

% Compute energies
E_main_lobe = sum(psf_patch(mask_circle).^2);
E_total_psf = sum(psf_patch(:).^2);
psf_ratio_dB = 10 * log10((E_total_psf-E_main_lobe) / E_total_psf); % simplified formula from lecture 3/4 p.5, slide 9

% Output result
fprintf('\n--- PSF Energy Analysis ---\n');
fprintf('Main lobe energy (within 2.5λ): %.2f\n', E_main_lobe);
fprintf('Total energy (10x10 mm region): %.2f\n', E_total_psf);
fprintf('Energy ratio (main lobe / total) in dB: %.2f dB\n', psf_ratio_dB);

%% Visualization
FontSz = 14;
set(groot, 'DefaultAxesFontSize', FontSz);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', FontSz);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', FontSz);       % Sets default font size for text elements

% Plot
figure('Position',[100, 100, 800, 800]);
imagesc(lateral_range(x_psf_range)*1000, depth_range(z_psf_range)*1000, 20*log10(psf_patch));
colormap(gray);
axis image;
colorbar;
title(sprintf('PSF Region at 30 mm (2.5λ = %.2f mm)', r_circle*1000));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
clim([-60, 0]);
hold on;
contour(lateral_range(x_psf_range)*1000, depth_range(z_psf_range)*1000, mask_circle, [1 1], 'Cyan', 'LineWidth', 1.5);

% Add contrast annotation
text(lateral_range(x_psf_range(1))*1000 + 1, depth_range(z_psf_range(1))*1000 + 1, ...
    sprintf('Main lobe ratio: %.2f dB', psf_ratio_dB), ...
    'Color', 'white', 'FontSize', 12, 'BackgroundColor', 'black', ...
    'EdgeColor', 'white', 'Margin', 5);

%% Show both PSF window and main lobe in full image

% Plot full image
figure('Position',[950, 100, 800, 800]);
imagesc(lateral_range*1000, depth_range*1000, 20*log10(rf_data));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
t = title('Contrast - Point Spread Function');
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.02;           % Move title higher (default is around 1.01)
axis ij equal tight;
clim([-60, 0]);
cb = colorbar; cb.Label.String = 'Received pressure (dB)';
cb.FontSize = FontSz;
colormap(gray);

% Draw 10x10 mm rectangle
hold on;
rectangle('Position', ...
    [lateral_range(xc_idx - half_win)*1000, ...
     depth_range(zc_idx - half_win)*1000, ...
     psf_window_size*1000, psf_window_size*1000], ...
    'EdgeColor', 'magenta', 'LineWidth', 2, 'LineStyle', '-');

% Overlay main lobe as a circle
theta = linspace(0, 2*pi, 200);
circle_x = lateral_psf*1000 + r_circle*1000 * cos(theta);
circle_z = depth_psf*1000 + r_circle*1000 * sin(theta);
plot(circle_x, circle_z, 'c', 'LineWidth', 1.5);

% Add legend text
text(lateral_range(1)*1000 + 5, 5, ...
    sprintf('Cystic Resolution: %.2f dB', psf_ratio_dB), ...
    'Color', 'white', 'FontSize', 12, 'BackgroundColor', 'black', ...
    'EdgeColor', 'white', 'Margin', 5);

