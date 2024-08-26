function [h] = waterDepth_f(h_toe,bathymetry)
% Description: Obtain an array of water depth along the cross-shore direction 
%   for a seawall fronted with a constant slope foreshore and a vegetated zone.

%   Input variables:
    %   1: h_toe = Water depth at the toe of seawall [m]. 
    %   2: bathymetry = from bathymetry_f.m. Only the second column is used. 
%   Output variables:
    %   1: h = array of water depth from offshore towards land [m]. 

% 1.0: Compute water depth array (h)
h = h_toe - bathymetry(:,2);

end