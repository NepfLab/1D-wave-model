function [K_amp_plotting] = drag_components_f(aw_list)

% Description: Given the wave progression matrix, compute the contribution 
%   of each of the 4 wave dissipation mechanisms.

%   Input variables:
    %   1: aw_list = array output from wave_model_f.m.
%   Output variables:
    %   1: K_amp_plotting = struct of vegetation characteristics.


%% 1.0: Initial set up to store results
aw_inc = aw_list(2:end,3) - aw_list(1:end-1,3);

K_v_vec = aw_list(1:end-1,6); 
K_sh_vec = aw_list(1:end-1,7); 
K_br_vec = aw_list(1:end-1,8); 
K_bed_vec = aw_list(1:end-1,9); 

% 1.1: shoaling coefficient may be greater than one
K_sh_inv_vec = NaN(size(K_sh_vec)); 
K_sh_binary = zeros(size(K_sh_vec)); 

for row = 1:length(K_sh_vec)
    if K_sh_vec(row) > 1
        K_sh_inv_vec(row) = 1 / K_sh_vec(row);
        K_sh_binary(row) = -1;
    else
        K_sh_inv_vec(row) = K_sh_vec(row);
        K_sh_binary(row) = 1;
    end
end


%% 2.0: Breakdown contributions from each of the four components
% 2.1: compute common demoninator
K_matrix = cat(2, K_v_vec, K_sh_inv_vec, K_br_vec, K_bed_vec); 
K_matrix_one_minus = 1 - K_matrix; % reduction factor

% 2.2: compute proportion of each K
sum_vec = sum(K_matrix_one_minus,2); 
K_prop = K_matrix_one_minus ./ sum_vec;

% 2.3: adjust for increasing amplitude from shoaling
K_prop(:,2) = K_prop(:,2).* K_sh_binary;

% 2.4: amplitude contribution by each K
sum_prop = sum(K_prop, 2);
x_vec = aw_inc ./ sum_prop; 
K_amp = K_prop .* x_vec;

% 2.5: cummulative amplitude reduction by each K
K_amp_cumm = cumsum(K_amp,1);

K_amp_plotting =  cat(1, zeros(1,4), K_amp_cumm);


end
