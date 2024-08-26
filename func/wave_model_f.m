function [aw_toe, aw_list, aw_table] = wave_model_f(h,vegArray,total_steps,dx,bool_sh,bool_v,aw_0,T)
% Description: 1-D wave model of wave propagation from offshore towards a amrsh fronted seawall.  

%   Input variables:
    %   1: h = array of water depth from offshore to land with same length as variable 'total_steps' [m].
    %   2: vegArray = array of veg_type at each step [m].
    %   3: total_steps = total number of spatial steps to model.
    %   4: dx = cross-shore spatial step between each iteration [m]. 
    %   5: bool_sh = boolean of 0 = no shoaling, 1 = with shoaling. 
    %   6: bool_v = boolean of 0 = no vegetation drag, 1 = with vegetation drag.
    %   7: aw_0 = offshore wave amplitude [m], assumed to be half of wave height. 
    %   8: T = % offshore wave period [s]
%   Output variables:
    %   1: aw_toe = wave amplitude at the toe of the seawall, which is the model final wave amplitude [m].
    %   2: aw_list = array form of wave modeling results for each iteration.
    %   3: aw_table = table form of wave modeling results for each iteration with column headers. 


%% 1.0: Intialize model parameters
dist_from_seawall = linspace(dx*total_steps,0,total_steps+1);

aw_i = aw_0; % initialise first iteration wave amplitude [m]
w = 2*pi()/T; % angular wave frequency [rad/s]
g = 9.81; % gravitational acceleration [ms^-2]
pho = 1000; % water density [kg/m3]

%%% 1.1: tabulate k, wavelength, c_g and wave state for each water depth
k_list = zeros(total_steps+1, 1);
wavelength = zeros(total_steps+1, 1);
c_g_list = zeros(total_steps+1, 1);
wave_state = zeros(total_steps+1, 1);

for i = 1:total_steps+1
    k_i = wave_number_f(T,h(i)); % compute wave number [m^-1] 
    k_list(i) = k_i; % add k value to k_list
    wavelength(i) = 2*pi()/k_i; % compute wavelengths [m] 
    c_g_list(i) = 1/2 * (1 + 2*k_i*h(i)/sinh(2*k_i*h(i))) ...
        * (g/k_i*tanh(k_i*h(i)))^(1/2); % compute group velocity [ms^-1]

    % 1.2: classify the wave state
    if h(i) < wavelength(i)/20
        wave_state(i) = 1; % shallow water wave
    elseif h(i) > wavelength(i)/2
        wave_state(i) = -1; % deep water wave
    else
        wave_state(i) = 0; % intermediate water wave
    end
end

% 1.3: find first U_w
if isempty(vegArray(1))
    z = h(1); % if the first position has no vegetation, so use water depth for z. 
else 
    % load vegetation params to find degree of submergence. 
    veg_params = vegParams_f(vegArray(1)); 
    l_s = veg_params.l_s;

    if h(1) >= l_s 
        z = l_s; % if marsh submerged, use full length
    else
        z = h(1); % if marsh submerged, use water depth
    end
end

k_i = k_list(1); 
h_i = h(1);

U_w_i = 2*pi()/T * aw_i * (cosh(k_i*z))/(sinh(k_i*h_i)); % wave orbital velocity [ms^-1] to be used in Cauchy number


%% 2.0 Compute wave amplitude modeling

aw_list = NaN(total_steps+1,11); % create list to store distances, water depth, aw, aw ratio, wave velocity, K_v, K_sh, state
aw_list(:,1) = dist_from_seawall; % x steps
aw_list(:,2) = h; % water depths
aw_list(1,3) = aw_0; 
aw_list(1,4) = NaN; % K_D_br, energy dissipation per broken wave
aw_list(1,5) = U_w_i; % wave orbital velocity
aw_list(1,6) = NaN; % K_v
aw_list(1,7) = NaN; % K_sh 
aw_list(1,8) = NaN; % K_br
aw_list(1,9) = NaN; % K_bed 
aw_list(:,10) = wave_state; % deep, intermediate, or shallow wave

for step = 1:total_steps
    A_w = U_w_i / w; % wave orbital excursion to be used for length ratios 
    k_i = k_list(step); 
    h_i = h(step);
    cg_i = c_g_list(step);
    wave_state_i = wave_state(step); 
    wavelength_i = wavelength(step);
    veg_type = vegArray(step);

    % 2.1: Skip step conditions
    if h(step+1) <= 0 % if there are no water. 
        disp('warning: No water depth at next step.') 
        aw_i_1 = NaN; % no value for wave height for no presence of water. 
        continue
    end

    %% 3.0: Compute K_v vegetation drag coef
    if veg_type == 'NoVeg' | bool_v == 0
        K_v = 1;
    else 
        %%% 3.1: Load vegetation parameters for the step
        veg_params = vegParams_f(veg_type); 
    
        lambda_p = veg_params.lambda_p;
        l_l = veg_params.l_l;
        l_s = veg_params.l_s;
        N_s = veg_params.N_s;
        C_s = veg_params.C_s;
        N_l = veg_params.N_l;
        b = veg_params.b;
        ElIl = veg_params.ElIl;
        EsIs = veg_params.EsIs;
        D = veg_params.D;
        Es = veg_params.Es;
        Is = veg_params.Is;
    
        %%% 3.2: Compute vegetation derived parameters
        % Compute length ratio
        L_l = l_l/A_w; % length ratio of leaf
        L_s = l_s/A_w; % length ratio of stem
        
        % Compute CD coef of drags
        KC_l = U_w_i*T/b; 
        KC_s = U_w_i*T/D; 
    
        if KC_l <= 10
            CD_l = 16*KC_l^(-0.52);
        else
            CD_l = max(1.95,10*KC_l^(-1/3));
        end
        
        if KC_s <= 11
            CD_s = 0.19*KC_s + 0.2;
        elseif KC_s <= 25
            CD_s = 7.6*KC_s;
        else 
            CD_s = max(1,2.9*KC_s^(-0.2));
        end
    
        K_l = 1; % for simplicity. K_l = 0.94 +- 0.06
        K_s = K_l*(CD_l/CD_s)^(1/4);
    
        % Compute cauchy numbers
        A_s = D * l_s; % frontal area of circular stem
        A_l = b * l_l; % frontal area of flat leaf
    
        Ca_s = pho*A_s*U_w_i^2 / (EsIs/l_s^2); % stem cauchy number
        Ca_l = pho*A_l*U_w_i^2 / (ElIl/l_l^2); % leaf cauchy number
    
        % Compute beta: Change in wave orbital velocity within the plant canopy
        if lambda_p <= 0.06
            beta = 1;
        elseif h_i - l_s < 0 % if marsh is submerged
            beta = 1; % check equation 18 in Lowe et. al. (2005) % NOT SURE
        else
            beta = 1/(1-lambda_p);
        end
    
        %%% 3.3: Compute K_v vegetation dissipation coef
    
        if l_l == 0
            K_D = 8/(9*pi())*k_i*beta^3 ...
                * ((sinh(k_i*l_s))^3 + 3*sinh(k_i*l_s)) / ((sinh(2*k_i*h_i)+2*k_i*h_i)*sinh(k_i*h_i))...
                * N_s * (...
                + CD_s*D*K_s*(Ca_s*L_s)^(-1/4)); % damping coefficient
        else
            K_D = 8/(9*pi())*k_i*beta^3 ...
                * ((sinh(k_i*l_s))^3 + 3*sinh(k_i*l_s)) / ((sinh(2*k_i*h_i)+2*k_i*h_i)*sinh(k_i*h_i))...
                * N_s * ((l_l/l_s*C_s*N_l*CD_l*b*K_l*(Ca_l*L_l)^(-1/4))...
                + CD_s*D*K_s*(Ca_s*L_s)^(-1/4)); % damping coefficient
        end
    
        K_v = 1 / (1+K_D*aw_i*dx); % vegetation dissipation coefficient
    end

    %% 4.0: Compute K_sh shoaling coef
    if bool_sh == 0
        K_sh = 0;
    elseif bool_sh == 1
        if isnan(c_g_list(step+1)) % no wave group velocity
            K_sh = 0 ; % no shoaling
        else
            K_sh = (c_g_list(step)/c_g_list(step+1))^(1/2); % shoaling coefficient
        end
    else
        disp("Invalid bool_sh input.")
    end
    
    
    %% 5.0: Compute K_breaking depth-induced wave breaking dissipation coef (Battjes and Janssen, 1978)

    %%% 5.1: Find breaker height (Miche's breaking criteria)
    if wave_state_i == -1 % deep water wave
        H_m = 0.142 * wavelength_i;
    elseif wave_state_i == 0 % intermediate water wave
        H_m = 1/7 * wavelength_i * tanh(2*pi()*h_i / wavelength_i);
    else % shallow water wave
        H_m = 0.88 * h_i;
    end

    %%% 5.2: Probability of wave breaking Q_b:
    if 2*aw_i >= H_m
        Q_b = 1;
    else
        Q_b_list = NaN(2,100001);
        Q_b_list(1,:) = linspace(0,1,size(Q_b_list,2));
        for index = 1:length(Q_b_list)
            Q_b_list(2,index) = abs(((1-Q_b_list(1,index)) / (log(Q_b_list(1,index))) + (2*aw_i/H_m)^2));
        end
        [M,I] = min(Q_b_list(2,:));

        Q_b = Q_b_list(1,I);
    end

    %%% 5.3: Ratio of dissipation
    K_D_br = 4 * aw_i / (T*h_i*cg_i);

    if K_D_br > 1
        K_D_br
        Q_b
        K_D_br = 1; % Upper bound for ratio of dissipation
        disp('Wave completely broken');
    end
    K_br = sqrt(1 - K_D_br * Q_b);

    %% 6.0: Compute K_bed bed friction coef (Chapter 9 P.268-269, Section 9.2.2 of Dean & Dalrymple, 1991)

    d_90 = 0.5; % sand diameter for which 90% sand is finner than d_90 [mm] 
    k_e = 2*d_90; 

    f = 0.1*(k_e/A_w)^(3/4); 
    if k_e/A_w <= 0.02
        disp("Invalid f assumption")
        disp(k_e/A_w)
    end

    K_D_bed = 2*f/(3*pi()) * k_i^2 / ((sinh(2*k_i*h_i)+2*k_i*h_i)*sinh(k_i*h_i)); 

    K_bed = 1/(1+K_D_bed*aw_i*dx); % bed friction dissipation coefficient


    %% 7.0: Compute new wave amplitude at next time step aw_i_1

    aw_ratio = K_v * K_sh * K_br * K_bed;
    aw_i_1 = aw_i * aw_ratio;
    

    %% 7.1: Save amplitude results
    aw_list(step,6) = K_v; % K_v
    aw_list(step,7) = K_sh; % K_sh
    aw_list(step,8) = K_br; % K_br
    aw_list(step,9) = K_bed; % K_bed
    aw_list(step+1,3) = aw_i_1; % Wave amplitude for next step. 
    aw_list(step,4) = K_D_br; % K_D_br, energy dissipation per broken wave

    aw_i = aw_i_1; % reset for next amplitude


    %%% 7.2: compute wave orbital velocity at next step i+1
    if vegArray(step+1) == "NoVeg"
        z = h(step+1); % no vegetation, so use water depth
    else
        % load vegetation params
        veg_params = vegParams_f(vegArray(step+1)); 
        l_s = veg_params.l_s;
    
        if h(step+1) >= l_s 
            z = l_s; % if marsh submerged, use full length
        else
            z = h(step+1); % if marsh submerged, use water depth
        end
    end

    k_i = k_list(step+1); 
    h_i = h(step+1);

    U_w_i = 2*pi()/T * aw_i * (cosh(k_i*z))/(sinh(k_i*h_i)); 
    
    aw_list(step+1,5) = U_w_i; % save wave orbital velocity at step i+1


end

%%% 8.0: Tabulate modeling output results
aw_list(:,11) = (aw_list(1,3)-aw_list(:,3))/aw_list(1,3)*100; % compute percent wave height reduction

aw_table = array2table(aw_list);

aw_table.Properties.VariableNames = {'Distance from aw_0 [m]','Water depth [m]', ...
           'Wave amplitude [m]','K_D_br', 'Wave velocity [m/s]','K_v','K_sh', ...
           'K_br', 'K_bed', 'Wave state', '"% aw reduction'};

aw_toe = aw_list(size(aw_list,1),3); % wave amplitude at the seawall


end