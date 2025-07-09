%% Init
addpath(genpath('.\code'))
addpath(genpath('.\data'))
addpath(genpath('.\Field_II_ver_3_30_windows'))
field_init();
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% displacement = lambda/64; % inter frame displacement of the tissue
displacement = lambda/16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scanning parameters
% Define how many transmissions
no_transm = 128;
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

%% Create the impulse response
impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(max(size(impulse_response)))';
excitation = sin(2*pi*f0*(0:1/fs:1/f0));
% Define the apodization vector
apo_vector = hanning(Active_elements);

%% Field II related initialization

% Create a linear array transmit aperture
emit_aperture = xdc_linear_array(N_elements, width, element_height, kerf, 1, 1, focus);

% Set the impulse response of the transmit aperture
xdc_impulse(emit_aperture, impulse_response);

% Set the excitation
xdc_excitation(emit_aperture, excitation);

% Create a linear array receive aperture
receive_aperture = xdc_linear_array(N_elements, width, element_height, kerf, 1, 1, focus);

% Set the impulse response of the receive aperture
xdc_impulse(receive_aperture, impulse_response);

% Set the sampling rate
set_field('fs',fs);

%% Create the phantom
% Cyst phantom
[phantom_positions, phantom_amplitudes] = mod_cyst_phantom(3/1000, [0, 30]/1000, 10000); % radius, pos[x,z], number of scatteres

% Make a point target phantom
% phantom_positions = [0 0 10; 0 0 20; 0 0 30; 0 0 40; 0 0 50]/1000;
% phantom_amplitudes = ones(5,1)*1e6;

%% Perform STA scanning
clear image_LRI times
% initialization
times = zeros(N_elements,1);
apo_vector_tx_matrix = zeros(N_elements,no_transm);
% iterate for number of transmissions
for i = 1:no_transm
    % Set transmit aperture apodization
    actual_apo_vector = zeros(N_elements,1);
    actual_apo_vector( (1:Active_elements) + i-1) = apo_vector;
    apo_vector_tx_matrix(:,i) = actual_apo_vector;
    
    % Set the focus and apodization for this direction
    xdc_center_focus(emit_aperture, [x, 0, 0]);
    xdc_focus(emit_aperture, t0, [x, 0, z_focus]);
    xdc_center_focus(receive_aperture, [x, 0, 0]);
    xdc_focus(receive_aperture, focus_time, [x, 0, 5]);
    xdc_apodization(emit_aperture, t0, actual_apo_vector');
    xdc_apodization(receive_aperture, t0, ones(N_elements,1)');
    
    % Calculate the received response using calc_scat_multi
    [v, t1] = calc_scat_multi(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    
    % Store the result with zero padding
    times(i) = t1;
    v2 = cat(1, zeros(round(t1*fs), N_elements), v);
    image_LRI(1:length(v2),:,i) = v2;

    % Move the beam
    x = x + d_x;

    % Tissue movement lateral
    % Point phantom, points are shifted to the right
    phantom_positions(:,1) = i*displacement;
    % Cycst phantom, all scatterers are shifted to the right
    phantom_positions(:,1) =  phantom_positions(:,1) + displacement;
    
    % Tissue movement axial
    % phantom_positions(:,3) = phantom_positions(:,3) + displacement;

end


%% Compute HRI by summing LRI / Normalize
image_HRI = sum(image_LRI, 3);
image_HRI = image_HRI / max(image_HRI, [], 'all');
image_LRI = image_LRI / max(image_LRI, [], 'all');
z = [0, size(image_LRI,1)/fs * c/2];

%% Visualize 
set(groot, 'DefaultAxesFontSize', 12);       % Sets default font size for axes
set(groot, 'DefaultColorbarFontSize', 12);   % Sets default font size for colorbar
set(groot, 'DefaultTextFontSize', 12);       % Sets default font size for text elements
cnt = 0;
% LRI
figure('Position', [100, 100, 1000, 400]);
for i = [1,64,128] % show 3 transmissions separately
    cnt = cnt+1;
    subplot(1,3,cnt) 
    imagesc([0, image_width*1000], [0, z*1000], 20*log10(abs(hilbert(image_LRI(:,:,i)))));
    xlabel('Lateral (mm)');
    ylabel('Depth (mm)');
    axis ij equal tight;
    clim([-60, 0]);
    cb = colorbar; cb.Label.String = 'Received pressure (dB)'; colormap(gray); hold on;
end
sgtitle('RF data of emission #1, #64, #128');
% HRI
figure;
imagesc([0, image_width*1000], [0, z*1000], 20*log10(abs(hilbert(image_HRI))));
xlabel('Lateral (mm)');
ylabel('Depth (mm)');
title('High Resolution Image');
axis ij equal tight;
clim([-60, 0]);
cb = colorbar; cb.Label.String = 'Received pressure (dB)'; colormap(gray);

