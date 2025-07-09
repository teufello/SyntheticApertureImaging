%% Init
addpath(genpath('.\code'))
addpath(genpath('.\data\simulation'))
addpath(genpath('.\Field_II_ver_3_30_windows'))
field_init();
clear all; close all;

%% Load acquired point RF data
load group3_point_noMotion.mat  % size 5952x192x128
% load group3_cyst_noMotion.mat

LRI = double(raw(91:end, :, :));
LRI_shifted = zeros(size(raw(91:end, :, :)), 'like', raw);

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 1480; % 1480 for point scatters, 1540 for the rest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = 7e6;
fs = 62.5e6;
lambda = c/f0; % wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displacement = lambda/64; % interframe displacement of the tissue 
displacement = lambda/16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % Spacing between virtual sources = 0.23/1000 (m)
% Define how many transmissions
no_transm = size(raw,3);
% Compute the image width
image_width = (no_transm-1)*pitch;
% Compute the position of the left most line
x = -image_width/2;

% Iterate over LRIs
for i = 1:size(raw, 3)
    shift_amount = (i - 1) * displacement / pitch; % Compute lateral shift in elements
    
    % Interpolate RF data to new lateral positions
    for j = 1:size(LRI, 1) % Loop over depth samples
        LRI_shifted(j, :, i) = interp1(1:size(LRI, 2), LRI(j, :, i), ...
                                       (1:size(LRI, 2)) - shift_amount, 'cubic', 0);
    end
end

% LRI/HRI
HRI = sum(LRI_shifted, 3);
HRI = HRI / max(HRI, [], 'all');
LRI = LRI_shifted;
z = [0, size(LRI, 1)/fs * c/2];

% Visu
set(groot, 'DefaultAxesFontSize', 12);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', 12);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', 12);       % Sets default font size for text elements
figure;
imagesc([0, image_width*1000], [0, z*1000], 20*log10(abs(hilbert(HRI))));
xlabel('Lateral (mm)', 'FontSize', 18);
ylabel('Depth (mm)', 'FontSize', 18);
title('High Resolution Image');
axis ij equal tight;
ax = gca;
ax.FontSize = 18;
clim([-60, 0]);
cb = colorbar('FontSize', 18); cb.Label.String = 'Received pressure (dB)'; colormap(gray);
