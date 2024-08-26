function [bathymetry] = bathymetry_f(offshore_dist,flat_width,dx,slope)
% Description: Obtain the bathymetry profile for a seawall fronted with a 
%   constant slope foreshore and a vegetated zone.

%   Input variables:
    %   1: offshore_dist = offshore distance measured from seawall [m] = 2000 
    %   2: flat_width = length of vegetation considered [m] = 5
    %   3: dx = incremental step size [m] 
    %   4: slope = ground slope of foreshore (dy/dx) 
%   Output variables:
    %   1: bathymetry = column 1 with distance from seawall and column 2 with 
    %       an array of bathymetry. 


% 1.0: Get number of steps for iterative calculation
total_steps = fix(offshore_dist/dx); 
flat_steps = fix(flat_width/dx); 
offshore_to_flat_steps = total_steps - flat_steps;

% 2.0: Get array of distance from seawall
dist_from_seawall = linspace(dx*total_steps,0,total_steps+1);
% dist_from_seawall = total_steps:-dx:0;

% 3.0: Get bathymetry z array from foreshore
z_0 = - slope * dx * offshore_to_flat_steps; % negative as land elevation becomes more negative further seawards

if flat_steps == 0 % if there is no flat surface infront of seawall.
    z_slope = linspace(z_0,0,offshore_to_flat_steps+1);
    z = z_slope;
else
    z_flat = zeros(1,flat_steps); % Set datum for bed elevation be 0 m at seawall.
    z_slope = linspace(z_0,z_flat(end),offshore_to_flat_steps+1);
    z = cat(2,z_slope,z_flat);
end

% 4.0: Concatenate 'dist_from_seawall' and 'z' into variable 'bathymetry'.
bathymetry = [dist_from_seawall', z'];

end