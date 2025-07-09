%% Init
clear all; close all;

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool('local', 3);  % max. 4
end

%% Define the parameters of the transducer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 100e6; % 62.5e6 for lab scanner, 100e6 for simulation
c = 1540; % 1480 for acq point scatters, 1540 for the rest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_elements = 192;
Active_elements = 64;
f0 = 7e6; % center frequency
lambda = c/f0; % wavelength
kerf = 0.03/1000;
width = 0.2/1000;
pitch = width + kerf; % Spacing between virtual sources = 0.23/1000 (m)
F_recv = 1;
F_transm = -1;
pixel_size = 0.1/1000; %  pixel size laterally/axially 

%%%%%%%%%%
%% Load %%
%%%%%%%%%%
% Simulation
% point phantom
% load image_LRI_point_phantom.mat
% moving phantom
% load image_LRI_motion_lambda16.mat
% load image_LRI_motion_lambda64.mat
% load image_LRI_sim_point_lambda32.mat
% load image_LRI_sim_point_lambda4.mat
load image_LRI_sim_point_lambda2.mat
rf_data  = image_LRI;

% Acquired data
% point phantom
% load group3_point_noMotion.mat
% rf_trunc = double(raw(91:end,:,:)); % truncate data due to acquisition

% moving phantom
% load image_LRI_acq_point_lambda64.mat
% load image_LRI_acq_point_lambda16.mat

% rf_trunc = double(LRI); % already truncated 

% Filter
f_low = 6e6;          % Lower cutoff frequency in Hz
f_high = 8e6;         % Upper cutoff frequency in Hz
filter_order = 128;   % Filter order 
% Normalize cutoff frequencies
wn = [f_low, f_high] / (fs / 2);
% Design FIR bandpass filter
fir_coeff = fir1(filter_order, wn, 'bandpass');
% rf_data = filtfilt(fir_coeff, 1, rf_trunc);
% % freqz(fir_coeff, 1, 1024, fs); % validate filter 

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
t0 = 0;
z_focus = (Active_elements*pitch)/F_transm;
% Define electronic focus
focus_time = 0;

% Define array
arrayPos = (-N_elements/2+0.5:N_elements/2-0.5)*pitch; % Make a pitched array
arrayPos = [arrayPos.', zeros(N_elements,1)]; % the z-coordinate of the array pos
% Define virtual sources
VS = [(-63.5:63.5)' * pitch, zeros(no_transm,1)]; % Using 0 for z-axis
% Apodization
apod = hanning(N_elements);

% Create the pixel map
num_depth = 70/1000/pixel_size;
num_lateral_upsampled = image_width/pixel_size*10; % Upsampling factor (10 delivers smooth results)
depth_range = linspace(0.1, 70, num_depth) / 1000; % Convert
lateral_range = linspace(-image_width/2, image_width/2, num_lateral_upsampled); 
[X,Z] = meshgrid(lateral_range, depth_range);
pixelMap = cat(3, X, Z);

%% Compute upsampled LRI for points
% Prepare scanning on point phantoms only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lines = int16([10.2, 20.2, 30.2, 40.3, 50.3]/1000/pixel_size); % [10.2, 20.2, 30.2, 40.3, 50.3]; [27.9, 52.4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRI = zeros(size(pixelMap,1),size(pixelMap,2),no_transm);
% Apply hilbert transform at once
rf_data_hilbert = hilbert(rf_data);

parfor i = 1:no_transm
    LRI(:,:,i) = dyn_imageFormation(rf_data_hilbert(:,:,i),pixelMap, lines, VS(i,:), arrayPos, c, fs);
    fprintf("Loop no.: %d of %d\n", i, no_transm);
end

HRI = sum(LRI,3);
HRI = abs(HRI) / max(abs(HRI), [], 'all');

%% Plot
set(groot, 'DefaultAxesFontSize', 12);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', 12);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', 12);       % Sets default font size for text elements

num_lines = numel(lines);
figure('Position', [100, 100, 800, 900]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgtitle('\bf Acquired point, tissue motion (\lambda/16)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:num_lines
    row_idx = lines(idx);
    profile = HRI(row_idx, :);
    profile_dB = 20*log10(profile / max(profile));  % dB scale

    % Find FWHM - considering only the main lobe
    half_max = -6; % dB
    [~, peak_idx] = max(profile_dB); % Find peak of main lobe
    
    % Find nearest left index below -6 dB
    left_idx = find(profile_dB(1:peak_idx) < half_max, 1, 'last');
    
    % Find nearest right index below -6 dB
    right_idx = find(profile_dB(peak_idx:end) < half_max, 1, 'first') + peak_idx - 1;
    
    % Sanity check for crossing existence
    if ~isempty(left_idx) && left_idx < peak_idx && ~isempty(right_idx) && right_idx > peak_idx
        % Interpolate precisely left position
        x_left = interp1(profile_dB(left_idx:left_idx+1), ...
                         lateral_range(left_idx:left_idx+1)*1000, half_max);
    
        % Interpolate precisely right position
        x_right = interp1(profile_dB(right_idx-1:right_idx), ...
                          lateral_range(right_idx-1:right_idx)*1000, half_max);
    
        % Calculate lateral resolution (FWHM in mm)
        lateral_res = abs(x_right - x_left);
        lambda_res = lateral_res/lambda/1000;
    else
        lateral_res = NaN; % invalid case, no crossings found
    end

    % Plot profile
    subplot(num_lines, 1, idx); 
    hold on; grid on;
    plot(lateral_range*1000, profile_dB, 'LineWidth', 1.5, 'DisplayName', 'Upsampled profile');
    yline(-6, '--r', 'LineWidth', 1, 'DisplayName', '-6 dB threshold');
    %xlim([-2 2]);
    ylim([-60 0]); yticks(-60:20:0);
    xlabel('Lateral (mm)');
    ylabel('Amplitude (dB)');
    title(sprintf('FWHM = %.2f mm (%.2f \\lambda) at %.1f mm depth', ...
    lateral_res, lambda_res, depth_range(row_idx)*1000));
    legend('Location', 'west');
end


%% Nested function
function LRI = dyn_imageFormation(rf_data, pixelMap, lines, VS, arrayPos, c, fs)

    si = size(pixelMap);
    LRI = zeros(si(1:2));
    t = (0:size(rf_data,1)-1) / fs;

    numElements = length(arrayPos);
    x = 1:numElements;
    [X, T] = meshgrid(x, t);
    X2 = repmat(x, [si(2), 1]);  % same x-axis for each row

    for i = lines
        delay_line = nan(si(2), numElements);  % [lateral points x elements]
        apod_line = zeros(si(2), numElements); % corresponding apodization

        for j = 1:si(2)
            pixel = [pixelMap(i,j,1), pixelMap(i,j,2)];
            [full_delay, apod_vec] = dyn_delayCal(pixel, [i,j], VS, arrayPos, c, si(2));

            delay_line(j,:) = full_delay';
            apod_line(j,:) = apod_vec';
        end

        % Interpolation with interp2: columns = channels, rows = time
        value = interp2(X, T, rf_data, X2, delay_line, 'cubic');

        % Apodization masks the relevant values
        value_weighted = value .* apod_line;

        % Row-wise summation → LRI
        LRI(i,:) = sum(value_weighted, 2, 'omitnan');
    end
end
