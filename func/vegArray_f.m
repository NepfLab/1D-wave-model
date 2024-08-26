function [vegArray] = vegArray_f(veg_types,bounds,bathymetry,bounds_metric)

% Description: Obtain an array of vegetation type for each computation step 
%   for a seawall fronted with a constant slope foreshore and a vegetated zone.

%   Input variables:
    %   1: veg_types = ix1 cell of string vectors for i vegetation types. 
    %   2: bounds = ix2 array of distance or elevation lower and upper bounds for i vegetation types
    %   3: bathymetry = array of bathymetry elevation at each iteration from offshore towards land [m].
    %   4: bounds_metric = "distance" or "elevation" to specify the units for bounds.
%   Output variables:
    %   1: vegArray = array of vegetation type at each incremental step.  


% 1.0: Compute vegetation array (vegArray)
vegArray = repmat("NoVeg", size(bathymetry,1), 1); % initialize with no vegetation.

if size(bounds,1) ~= 0 % check if there are vegetation. 
    for i = 1:size(bounds,1) % inclusive lower bound, exclusive upperbound.
        lowerBound = bounds(i,1);
        upperBound = bounds(i,2);
        if bounds_metric == "distance" % specify bounds by distance noted in first column of bathymetry. 
            vegArray(bathymetry(:,1) >= lowerBound & bathymetry(:,1) < upperBound) = veg_types(i);
        elseif bounds_metric == "elevation" % specify bounds by ground elevation. 
            vegArray(bathymetry(:,2) >= lowerBound & bathymetry(:,2) < upperBound) = veg_types(i);
        end

        %%%%% Note on 'bounds':
        %%%%% 1) If the first column of 'bathymetry' is descending, such as 
        %%%%%       distance from seawall, the bounds = [lowerBound upperBound] + dx.   
        %%%%% 2) If the first column of 'bathymetry' is ascending, such as 
        %%%%%       numbering from distance from outer most marsh edge following 
        %%%%%       wave propagation, the bounds = [lowerBound upperBound]. 
        
        %%%%% General note:
        %%%%% vegArray outputs the vegetation type for each wave modeling
        %%%%%       iteration. 
    end
end

end