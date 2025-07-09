%% Init
addpath(genpath('.\code'))
addpath(genpath('.\data'))
addpath(genpath('.\Field_II_ver_3_30_windows'))
clear all; close all;

%%%%%%%%%%
%% Load %%
%%%%%%%%%%
% load image_HRI_beamforming_sim_cyst.mat
% load image_HRI_beamforming_sim_cyst_lambda64.mat
% load image_HRI_beamforming_sim_cyst_lambda16.mat

load image_HRI_beamforming_acq_cyst.mat
% load image_HRI_beamforming_acq_cyst_lambda64.mat
% load image_HRI_beamforming_acq_cyst_lambda16.mat
rf_data  = abs(HRI);

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define cyst centers [x_mm, z_mm]
cyst_centers_mm = [ -0.8, 16.4; -0.2, 44.9];  % [ 0, 17; 0, 47] for acq cyst; [0, 30] for sim
no_transm = 128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf;
pixel_size = 0.1/1000;

image_width = (no_transm-1)*pitch;
x = -image_width/2;
d_x = pitch;

num_depth = 70/1000/pixel_size;
num_lateral = image_width/pixel_size;
depth_range = linspace(0.1, 70, num_depth) / 1000;
lateral_range = linspace(-image_width/2, image_width/2, num_lateral);

%%%%%%%%%%%%%%
%% Contrast %%
%%%%%%%%%%%%%%
% CNR
mm_per_pixel = pixel_size * 1000;
window_size = 2.5/1000;
window_size_pixel = round(window_size / pixel_size);
half_win = floor(window_size_pixel/2);

%% Visu
FontSz = 14;
set(groot, 'DefaultAxesFontSize', FontSz);
set(groot, 'DefaultColorbarFontSize', FontSz);
set(groot, 'DefaultTextFontSize', FontSz);
% Plot
figure('Position',[100, 100, 800, 800]);
imagesc(lateral_range*1000, depth_range*1000, 20*log10(rf_data));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tissue motion ({\lambda}/16)
t = title('Acquired cyst, no motion');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.02;           % Move title higher (default is around 1.01)
axis ij equal tight;
clim([-60, 0]);
cb = colorbar;
cb.Label.String = 'Received pressure (dB)';
cb.FontSize = FontSz;
colormap(gray);
hold on;

for i = 1:size(cyst_centers_mm,1)
    lateral_center = cyst_centers_mm(i,1) / 1000;
    depth_center = cyst_centers_mm(i,2) / 1000;

    [~, z_idx] = min(abs(depth_range - depth_center));
    [~, x_idx] = min(abs(lateral_range - lateral_center));

    z_range = (z_idx-half_win):(z_idx+half_win);
    x_range = (x_idx-half_win):(x_idx+half_win);

    window_cyst = rf_data(z_range, x_range);

    lateral_tissue = lateral_center + 5/1000;
    [~, xg_idx] = min(abs(lateral_range - lateral_tissue));
    xg_range = (xg_idx-half_win):(xg_idx+half_win);

    window_tissue = rf_data(z_range, xg_range);

    Mu_cyst = mean(window_cyst, 'all');
    Mu_tissue = mean(window_tissue, 'all');
    Sigma_cyst = std(window_cyst, 0, 'all');
    Sigma_tissue = std(window_tissue, 0, 'all');

    cnr = 20*log10(abs((Mu_cyst - Mu_tissue)) / sqrt(Sigma_cyst^2 + Sigma_tissue^2));
    fprintf('  CNR: %.2f dB\n\n', cnr);

    % Draw cyst window and label
    cyst_x = lateral_range(x_idx-half_win)*1000;
    cyst_y = depth_range(z_idx-half_win)*1000;
    rectangle('Position', [cyst_x, cyst_y, window_size*1000, window_size*1000], 'EdgeColor', 'cyan', 'LineWidth', 2);

    % Draw tissue window and label
    tissue_x = lateral_range(xg_idx-half_win)*1000;
    tissue_y = depth_range(z_idx-half_win)*1000;
    rectangle('Position', [tissue_x, tissue_y, window_size*1000, window_size*1000], 'EdgeColor', 'magenta', 'LineWidth', 2);
    
    % Add text box with contrast info
    text(-image_width*1000/2 + 8, cyst_centers_mm(i,2)+8, sprintf('CNR = %.2f dB', cnr), ...
         'Color', 'white', 'FontSize', FontSz, 'BackgroundColor', 'black', ...
         'Margin', 5, 'EdgeColor', 'white', 'HorizontalAlignment', 'left');
    hold on;
end

text(-image_width*1000/2 + 2.5, 3, '\color{cyan}⬜\color{white} Cyst ROI   \color{magenta}⬜\color{white} Tissue ROI', ...
     'FontSize', FontSz, 'BackgroundColor', 'black', 'EdgeColor', 'white', 'Margin', 5);