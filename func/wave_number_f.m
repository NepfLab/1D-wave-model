function [wave_number] = wave_number_f(wave_period,water_depth)
% Description: Iterative process to obtain wave number. 

%   Input variables:
    %   1: wave_period = Wave period (T) [s]
    %   2: water_depth = Water_depth (h) [m] 
%   Output variables: 
    %   1: wave_number = wave number (k) [m^-1]. 

%% 1.0: Initial set up before iterative solution. 
if water_depth <= 0
    disp('wave_number_f error: Negative water_depth.')
end

g = 9.81; % gravitational acceleration [ms^-2]

% initialise first guess for wave number with deep water assumption. 
rad_freq = 2*pi()/wave_period; % fixed
wave_number = rad_freq^2/g; % initialise first guess for wave number. 

%% 2.0: Iterate to compute wave number
error = 1; % store the difference between each wave number iteration. 
threshold = 0.000001; % threshold to end iteration on wave number. 

while(error)>threshold
    wave_number_iter = rad_freq^2 / (g*tanh(wave_number*water_depth)); % dispersion relationship. 
    error = abs(wave_number - wave_number_iter); 
    wave_number = wave_number_iter; 
end