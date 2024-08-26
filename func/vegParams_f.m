function [veg_params] = vegParams_f(veg_type)
% Description: Obtain the vegetation parameters of the specified veg_type as a structure data. 

%   Input variables:
    %   1: veg_type = string vegetation species name.
%   Output variables:
    %   1: veg_params = struct of vegetation characteristics.


% 1.0: Create a structure with the given parameters
    veg_params = struct();
    veg_params.lambda_p = NaN; % solid volume fraction
    veg_params.l_l = NaN; % length of leaf [m]
    veg_params.l_s = NaN; % height of marsh [m]
    veg_params.N_s = NaN; % number of stem per area [m^-2]
    veg_params.C_s = NaN; % sheltering coefficient
    veg_params.N_l = NaN; % number of leaves [stem^-1]
    veg_params.b = NaN; % leaf width; width at leaf base [m]
    veg_params.d = NaN; % leaf thickness [m]
    veg_params.D = NaN; % stem diameter [m]
    veg_params.El = NaN; % leaf Young's modulus [GPa]
    veg_params.Es = NaN; % stem Young's modulus [GPa]
    veg_params.Il = NaN; % leaf bending moment of inertia [m^4]
    veg_params.Is = NaN; % stem bending moment of inertia [m^4]
    veg_params.ElIl = NaN; % leaf Young's modulus * leaf bending moment of inertia [Nm^-2]
    veg_params.EsIs = NaN; % stem Young's modulus * stem bending moment of inertia [Nm^-2]

    if veg_type == "NoVeg"
        % If veg_type = NoVeg, skip processing
        return;
    end

    switch veg_type
        % 2.0: Standard vegetation species. 
        case 'Phragmites australis' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.3; 
            veg_params.l_s = 2; 
            veg_params.N_l = 13; 
            veg_params.N_s = 340; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.02; 
            veg_params.d = 0.000285; 
            veg_params.D = 0.009; 
            veg_params.El = 3; 
            veg_params.Es = 7.5; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Scirpus maritimus' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.3; 
            veg_params.l_s = 0.85; 
            veg_params.N_l = 6; 
            veg_params.N_s = 550;
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.007; 
            veg_params.d = 0.0003; 
            veg_params.D = 0.006; 
            veg_params.El = 3; 
            veg_params.Es = 0.8; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Scirpus tabernaemontani' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0; 
            veg_params.l_s = 1.5; 
            veg_params.N_l = 0; 
            veg_params.N_s = 550; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0; 
            veg_params.d = 0; 
            veg_params.D = 0.012;      
            veg_params.El = 3; 
            veg_params.Es = 0.125; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Spartina alterniflora' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.35; 
            veg_params.l_s = 1.3; 
            veg_params.N_l = 9; 
            veg_params.N_s = 295; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.011; 
            veg_params.d = 0.00046; 
            veg_params.D = 0.007;  
            veg_params.El = 3; 
            veg_params.Es = 0.74;             
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Spartina anglica'  % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.3; 
            veg_params.l_s = 0.65; 
            veg_params.N_l = 8; 
            veg_params.N_s = 655; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.011; 
            veg_params.d = 0.000335; 
            veg_params.D = 0.005;    
            veg_params.El = 3; 
            veg_params.Es = 0.4; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
    
        case 'Scirpus mariqueter' % % Zhang et al., 2021 P1 scirpus mariqueter 
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.38; 
            veg_params.N_l = 1; 
            veg_params.N_s = 2352; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.006; 
            veg_params.d = 0.0003; 
            veg_params.D = 0.0022; 
            veg_params.El = 0.06; 
            veg_params.Es = 0.12; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;            
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3.0: Vegetation species in dormant.
        case 'Phragmites australis (dormant)' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.5; 
            veg_params.N_l = 7; 
            veg_params.N_s = 80; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.01; 
            veg_params.d = 0.00024; 
            veg_params.D = 0.003; 
            veg_params.El = 3; 
            veg_params.Es = 5; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Scirpus maritimus (dormant)' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.5; 
            veg_params.N_l = 4; 
            veg_params.N_s = 100; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.002; 
            veg_params.d = 0.0002; 
            veg_params.D = 0.003; 
            veg_params.El = 3; 
            veg_params.Es = 0.4;             
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Scirpus tabernaemontani (dormant)' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0; 
            veg_params.l_s = 1; 
            veg_params.N_l = 0; 
            veg_params.N_s = 300; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0; 
            veg_params.d = 0; 
            veg_params.D = 0.009;   
            veg_params.El = 3; 
            veg_params.Es = 0.12;             
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Spartina alterniflora (dormant)' % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.3; 
            veg_params.N_l = 5;             
            veg_params.N_s = 170; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.005; 
            veg_params.d = 0.00025; 
            veg_params.D = 0.003; 
            veg_params.El = 3; 
            veg_params.Es = 0.08; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;            
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Spartina anglica (dormant)'  % Zhang et al., 2021 supplementary
            veg_params.lambda_p = 0.006; % solid volume fraction
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.3; 
            veg_params.N_l = 5; 
            veg_params.N_s = 430; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.006; 
            veg_params.d = 0.0002; 
            veg_params.D = 0.002;    
            veg_params.El = 3; 
            veg_params.Es = 0.1; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;            
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9);
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4.0: Vegetation species for Fig. 1 Bath and Hellegat data validation.
        % 4.1: Bath vegetation characteristics. 
        case 'Bath' % Generic (Scirpus maritimus) but lower bound
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.5; 
            veg_params.N_l = 4; 
            veg_params.N_s = 100; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.002; 
            veg_params.d = 0.0002; 
            veg_params.D = 0.003; 
            veg_params.El = 3; 
            veg_params.Es = 0.4; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            % veg_params.ElIl = 4.0*10^(-6); 
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9);            

        case 'Bath S1-S2'
            veg_params = vegParams_f('Bath');
    
            veg_params.l_s = 0.17; 
            veg_params.N_s = 144; 
            veg_params.D = 0.0087; 
            veg_params.lambda_p = 1.5*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
    
        case 'Bath S2-S3'
            veg_params = vegParams_f('Bath');
    
            veg_params.l_s = 0.15; 
            veg_params.N_s = 372;
            veg_params.D = 0.008; 
            veg_params.lambda_p = 2.8*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
        
        case 'Bath S3-S4'
            veg_params = vegParams_f('Bath');
    
            veg_params.l_s = 0.35; 
            veg_params.N_s = 1072; 
            veg_params.D = 0.0049; 
            veg_params.lambda_p = 7.1*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        % 4.2: Hellegat vegetation characteristics. 
        case 'Hellegat' % Generic (Spartina anglica, Zhang 2021 Supplementary lower bound)
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.3; 
            veg_params.N_l = 5; 
            veg_params.N_s = 430; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.006; 
            veg_params.d = 0.0002; 
            veg_params.D = 0.002; 
            veg_params.El = 3; 
            veg_params.Es = 0.1; 
            veg_params.Il = veg_params.b * veg_params.d^3/12;             
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9);             
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Hellegat S1-S2'
            veg_params = vegParams_f('Hellegat');
    
            veg_params.l_s = 0.2; 
            veg_params.N_s = 944; 
            veg_params.D = 0.003; 
            veg_params.lambda_p = 1.3*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
        
        case 'Hellegat S2-S3'
            veg_params = vegParams_f('Hellegat');
    
            veg_params.l_s = 0.29; 
            veg_params.N_s = 1136; 
            veg_params.D = 0.0034; 
            veg_params.lambda_p = 3*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
    
        case 'Hellegat S3-S4'
            veg_params = vegParams_f('Hellegat');
    
            veg_params.l_s = 0.27; 
            veg_params.N_s = 1520; 
            veg_params.D = 0.0027; 
            veg_params.lambda_p = 2.3*10^(-3); 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 


        % 5.0: Vegetation species for Fig. 2a sensitivity analysis.
        % 5.1: Bath vegetation characteristics. 
        case 'Bath S1-S2 low'
            veg_params = vegParams_f('Bath S1-S2');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Bath S2-S3 low'
            veg_params = vegParams_f('Bath S2-S3');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Bath S3-S4 low'
            veg_params = vegParams_f('Bath S3-S4');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Bath S1-S2 mid'
            veg_params = vegParams_f('Bath S1-S2');
            veg_params.N_l = 2; % number of leaves [stem^-1]

        case 'Bath S2-S3 mid'
            veg_params = vegParams_f('Bath S2-S3');
            veg_params.N_l = 2; % number of leaves [stem^-1]

        case 'Bath S3-S4 mid'
            veg_params = vegParams_f('Bath S3-S4');
            veg_params.N_l = 2; % number of leaves [stem^-1]

        % 5.2: Hellegat vegetation characteristics. 
        case 'Hellegat S1-S2 low'
            veg_params = vegParams_f('Hellegat S1-S2');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Hellegat S2-S3 low'
            veg_params = vegParams_f('Hellegat S2-S3');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Hellegat S3-S4 low'
            veg_params = vegParams_f('Hellegat S3-S4');
            veg_params.N_l = 0; % number of leaves [stem^-1]

        case 'Hellegat S1-S2 mid'
            veg_params = vegParams_f('Hellegat S1-S2');
            veg_params.N_l = 2; % number of leaves [stem^-1]

        case 'Hellegat S2-S3 mid'
            veg_params = vegParams_f('Hellegat S2-S3');
            veg_params.N_l = 2; % number of leaves [stem^-1]

        case 'Hellegat S3-S4 mid'
            veg_params = vegParams_f('Hellegat S3-S4');
            veg_params.N_l = 2; % number of leaves [stem^-1]


        % 6.0: Juniper Cove vegetation characteristics 
        case 'Juniper: Spartina alterniflora healthy'
            veg_params = vegParams_f('Spartina alterniflora');

        case 'Juniper: Spartina alterniflora dormant'
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.2; 
            veg_params.l_s = 0.3; 
            veg_params.N_l = 5;             
            veg_params.N_s = 170; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.005; 
            veg_params.d = 0.00025; 
            veg_params.D = 0.003;
            veg_params.El = 3; 
            veg_params.Es = 0.08; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;            
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9);

        case 'Juniper: Spartina patens healthy'
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.44; 
            veg_params.l_s = 0.76; 
            veg_params.N_l = 5.4; 
            veg_params.N_s = 752; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.0043; 
            veg_params.d = 0.00046; 
            veg_params.D = 0.0016;  
            veg_params.El = 3; 
            veg_params.Es = 1.5; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        case 'Juniper: Spartina patens dormant'
            veg_params.lambda_p = 0.006; 
            veg_params.l_l = 0.25; 
            veg_params.l_s = 0.34; 
            veg_params.N_l = 4; 
            veg_params.N_s = 563; 
            veg_params.C_s = 0.4; % sheltering coefficient -> Spartina alterniflora (Zhang, Nepf 2021)
            veg_params.b = 0.003; 
            veg_params.d = 0.00046; 
            veg_params.D = 0.0011;  
            veg_params.El = 3; 
            veg_params.Es = 0.13; 
            veg_params.Il = veg_params.b * veg_params.d^3 / 12; 
            veg_params.Is = pi() * veg_params.D^4 / 64;
            veg_params.ElIl = veg_params.El * veg_params.Il * 10^(9); 
            veg_params.EsIs = veg_params.Es * veg_params.Is * 10^(9); 

        otherwise
            disp('Invalid veg input')
    end

end