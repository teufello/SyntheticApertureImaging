% Create a computer model of a cyst phantom. The phantom contains
% five point targets separated by 5 mm and a 10 mm water filled cyst.
% All scatterers are situated in a box of (x,y,z)=(40,10,50) mm.
%
% Calling: [positions, amp] = cyst_phantom (N);
%
% Parameters: N - Number of scatterers in the phantom
%
% Output: positions - Positions of the scatterers.
% amp - amplitude of the scatterers.
%
% Version 1.1, March 22, 2011 by Joergen Arendt Jensen
% 
% Modified by:
% Andrea Matamoros
% Louis Teufel
% r: radius of cyst in mm
% pos: position of cyst [x,z] in mm

function [positions, amp] = mod_cyst_phantom (r,pos,N)
x_size = 15/1000; % Width of phantom [m]
y_size = 0/1000; % Transverse width of phantom [m]
z_size = 15/1000; % Height of phantom [m]
z_start = pos(2) - (z_size - r)/2; % Start of phantom surface [m];
% Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
% Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);
% Make the cyst and set the amplitudes to zero inside
inside = ( ((x-pos(1)).^2 + (z-pos(2)).^2) < r^2);
amp = amp .* (1-inside);

% Return the variables
positions=[x y z];
end