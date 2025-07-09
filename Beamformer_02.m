%% Init
addpath(genpath('.\code'))
addpath(genpath('.\data'))
addpath(genpath('.\Field_II_ver_3_30_windows'))
field_init();
clear all; close all;

%% Run in cluster
% delete(gcp('nocreate'))
% clust = parcluster("LouisTeufel_247239"); % using nodes
% if isempty(gcp('nocreate'))
%     p = parpool('local',20); % max. 20
% end

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local', 4);  % or whatever number of workers you want
end

%% Define the parameters of the transducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_data = true; % choose true: simulated data; false: acquired data
c = 1540; % 1480 for point scatterers, 1540 for the rest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim_data
    fs = 100e6; % 100e6 for simulation
else
    fs = 62.5e6; % 62.5e6 for lab scanner 
end
N_elements = 192;
Active_elements = 64;
f0 = 7e6; % Center freq in Hz
lambda = c/f0; % wavelength
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % Spacing between virtual sources = 0.23/1000 (m)
F_transm = -1;
pixel_size = 0.1/1000; %  pixel size laterally/axially 

%%%%%%%%%%
%% Load %% 
%%%%%%%%%%
%% Truncate and filter acquired data
if ~sim_data
% load group3_point_noMotion.mat
% load group3_cyst_noMotion.mat
% rf_trunc = double(raw(91:end,:,:)); % truncate data due to acquisition 

load image_LRI_acq_cyst_lambda64.mat % already truncated
% load image_LRI_acq_cyst_lambda16.mat % already truncated

% load image_LRI_acq_point_lambda64.mat % already truncated
% load image_LRI_acq_point_lambda16.mat % already truncated
% load image_LRI_0.mat

rf_trunc = double(LRI); % int16 to double

% Filter
f_low = 6e6;          % Lower cutoff frequency in Hz
f_high = 8e6;         % Upper cutoff frequency in Hz
filter_order = 128;   % Filter order 
% Normalize cutoff frequencies
wn = [f_low, f_high] / (fs / 2);
% Design FIR bandpass filter
fir_coeff = fir1(filter_order, wn, 'bandpass');
rf_data = filtfilt(fir_coeff, 1, rf_trunc);
% freqz(fir_coeff, 1, 1024, fs); % validate filter 
end

%% Load phantom
if sim_data
% Point phantoms
% load image_LRI_point_phantom.mat
% rf_data = image_LRI;

% Moving
load image_LRI_motion_lambda16.mat
% load image_LRI_motion_lambda64.mat
rf_data  = image_LRI;

% Cyst phantom 
SNR_dB = 20; % add noise
% load image_LRI_cyst_phantom.mat
% rf_data  = awgn(image_LRI, SNR_dB, 'measured');

% Moving
% load image_LRI_cyst_lambda16.mat
% load image_LRI_cyst_lambda64.mat
% rf_data  = awgn(image_LRI, SNR_dB, 'measured');
end

%% Parameters
% Define how many transmissions
no_transm = size(rf_data,3);
% Compute the image width
image_width = (no_transm-1)*pitch;
% Compute the position of the left most line
x = -image_width/2;
% Compute how much to display per line
d_x = pitch;
% Define transmit focus
z_focus = (Active_elements*pitch)/F_transm;

% Define array
arrayPos = (-N_elements/2+0.5:N_elements/2-0.5)*pitch; % Make a pitched array
arrayPos = [arrayPos.', zeros(N_elements,1)]; % the z-coordinate of the array pos
% Define virtual sources
VS = [(-63.5:63.5)' * pitch, ones(no_transm,1)*z_focus]; 
% Apodization
apod = hanning(N_elements);

% Create the pixel map
depth_mm = 70; % mm
lateral_mm = 30; % mm
num_depth = depth_mm/1000/pixel_size;
num_lateral = lateral_mm/1000/pixel_size;
% num_lateral = image_width/pixel_size;
depth_range = linspace(0.1, depth_mm, num_depth) / 1000; % Convert
lateral_range = linspace(-image_width/2, image_width/2, num_lateral); 
[X,Z] = meshgrid(lateral_range, depth_range);
pixelMap = cat(3, X, Z);

%%%%%%%%%%%%%%%%%
%% Beamforming %%
%%%%%%%%%%%%%%%%%
% Compute LRIs
LRI = zeros(size(pixelMap,1),size(pixelMap,2),no_transm);
rf_data_hilbert = hilbert(rf_data);
parfor i = 1:no_transm
    LRI(:,:,i) = dyn_imageFormation(rf_data_hilbert(:,:,i),pixelMap, VS(i,:), arrayPos, c, fs);
    fprintf("Loop no.: %d of %d\n", i, no_transm);
end
% Create HRI
HRI = sum(LRI,3);
HRI = abs(HRI) / max(abs(HRI), [], 'all');

%% Visu
set(groot, 'DefaultAxesFontSize', 14);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', 14);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', 14);       % Sets default font size for text elements
figure('Position',[100,100,600,800]);

% Plot
imagesc([-lateral_mm/1000/2*1000, lateral_mm/1000/2*1000], [0, depth_range*1000], 20*log10(HRI));
xlabel('Lateral (mm)', 'FontSize',14);
ylabel('Depth (mm)', 'FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = title('Acquired cyst, no motion', 'Interpreter', 'tex','FontSize',14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t.Units = 'normalized';         % Use normalized coordinates (0 to 1)
t.Position(2) = 1.02;           % Move title higher (default is around 1.01)
axis ij equal tight;
% ylim([20, 40]); % used for cyst phantom
clim([-60, 0]);
cb = colorbar; cb.Label.String = 'Received pressure (dB)'; cb.FontSize = 14; colormap(gray);
