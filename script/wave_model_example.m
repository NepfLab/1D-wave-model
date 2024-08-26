%% Script for 1-D wave model example

% Author(s): Ernie I. H. Lee, and Heidi Nepf, 2024.
% Email contact: eihlee@mit.edu. 

% Description: 1-D wave model example of wave propagation from offshore
%    towards a marsh fronted seawall. 

% Functions: All the custom MATLAB functions needed are stored in the 'func' folder. 
    % 1: bathymetry_f: Compute bathymetry array for the hypothetical marsh fronted seawall configuration.
    % 2: waterDepth_f: Compute water depth array.
    % 3: vegArray_f: Compute vegetation array.
    % 4: wave_model_f: Iterative 1-D wave modeling of wave propagation in the cross-shore direction. 
        % 4.1: wave_number_f: Compute wave number.
        % 4.2: vegParams_f: Retrieve vegetation parameters for the specific marsh species. 
    % 5: drag_components_f: Compute the contribution of each of the 4 wave dissipation mechanisms

%% 1.0 Initial set up
clear(); clc(); close all force;


%% 1.1 Set up folder directories
% Get the current directory
currentDir = pwd;

% Find the parent directory
[parentDir, ~, ~] = fileparts(currentDir);

% Specify the name of the new subfolder 'func', contains functions. 
func_folder = 'func';

% Add func subfolder containing functions into the file path. 
func_folder_Path = fullfile(parentDir, func_folder);

addpath(func_folder_Path); 

%% 2.0: 1-D Wave model

% 2.1: Model settings 
dx = 1; % incremental step size [m]
bool_sh = 1; % 0 and 1 for neglecting and including shoaling
bool_v = 1; % 0 and 1 for neglecting and including vegetation

%%% 2.2: seawall settings
offshore_dist = 100; % [m]
flat_width = 0; % width of flat platform in front of the seawall [m]
slope = 0.01; % ground slope for water depth []

%%% 2.3: hydrodynamic settings
h_toe = 3; % water depth at the seawall toe, or on the landward most point [m]
aw_0 = 1.5/2 ; % initialise first wave amplitude assuming wave amplitude is half of wave height [m]
T = 5; % wave period [s]

%%% 2.4: vegetation settings (choose either 3.3a or 3.3b) 
% %%% 2.4a: vegetation settings for 1 vegetation type
% vegTypes = {'Spartina Alterniflora'};
% veg_width = 50; % assume vegetation width is the same as flat_width [m]
% bounds = [0+dx,veg_width+dx];

%%% 2.4b: vegetation settings for more than 1 vegetation type
veg_types = {'Phragmites australis','Spartina alterniflora'};
bounds = [0,15; 15,30] + dx; % distance [m] from seawall for presence of vegetation by type. 
bounds_metric = "distance"; % define bounds by the iteration number. 

%%%%% Note on 'bounds':
%%%%% 1) If the first column of 'bathymetry' is descending, such as 
%%%%%       distance from seawall, the bounds = [lowerBound upperBound] + dx.   
%%%%% 2) If the first column of 'bathymetry' is ascending, such as 
%%%%%       numbering from distance from outer most marsh edge following 
%%%%%       wave propagation, the bounds = [lowerBound upperBound]. 

%%%%% General note:
%%%%% vegArray outputs the vegetation type for each wave modeling
%%%%%       iteration. 


%%% 2.5: Get bathymetry, water depth, veg array, and total_steps
bathymetry = bathymetry_f(offshore_dist,flat_width,dx,slope);
h = waterDepth_f(h_toe,bathymetry);
vegArray = vegArray_f(veg_types,bounds,bathymetry,bounds_metric);
total_steps = fix(offshore_dist/dx); % number of steps for iterative calculation


%%% 2.6: Predict wave amplitude set up
% with vegetation drag
bool_v = 1; % 0 and 1 for neglecting and including vegetation
[aw_toe, aw_list, aw_table] = wave_model_f(h,vegArray,total_steps,dx,bool_sh,bool_v,aw_0,T); 
aw_list_veg = aw_list; 

% without vegetation drag
bool_v = 0; % 0 and 1 for neglecting and including vegetation
[aw_toe, aw_list, aw_table] = wave_model_f(h,vegArray,total_steps,dx,bool_sh,bool_v,aw_0,T); 
aw_list_noveg = aw_list; 

%%% 2.7: Save aw_list to csv
aw_list_veg_csv = "aw_list_veg.csv";
aw_list_noveg_csv = "aw_list_noveg.csv";

writematrix(aw_list_veg,aw_list_veg_csv)
writematrix(aw_list_noveg,aw_list_noveg_csv)


%% 3.0 Visualize results

%%% 3.1 Load data
aw_list_veg_csv = "aw_list_veg.csv";
aw_list_noveg_csv = "aw_list_noveg.csv";

aw_list_veg = readmatrix(aw_list_veg_csv);
aw_list_noveg = readmatrix(aw_list_noveg_csv);

%%% 3.2 Create a new figure
figure

%%% 3.3: Left y-axis: Wave height plots
yyaxis left
p1 = area([bounds(1,1)-dx bounds(1,2)-dx], [110 110], 'FaceColor', 'g', 'FaceAlpha',0.1, 'DisplayName', veg_types{1});
hold on
p2 = area([bounds(2,1)-dx bounds(2,2)-dx], [110 110], 'FaceColor', 'g', 'FaceAlpha',0.3, 'DisplayName', veg_types{2});

p3 = plot(aw_list_veg(:,1), aw_list_veg(:,3)*2, '-', 'LineWidth', 3, ... 
        'DisplayName', 'With Vegetation'); % wave height

p4 = plot(aw_list_veg(:,1), aw_list_noveg(:,3)*2, '--', 'LineWidth', 3, ...
    'DisplayName', 'Without Vegetation'); % wave height

ylim([0 max(aw_list_noveg(:,3))*2+0.1]); % wave height
xlabel('Distance from seawall (m)', 'fontname', 'arial'); % wave height
ylabel('Predicted wave height (m)', 'fontname', 'arial'); % wave height

set(gca,'FontSize', 20, 'color','none');

%%% 3.4: Right y-axis: Bathymetry
yyaxis right
p5 = plot(aw_list_veg(:,1),bathymetry(:,2),'LineWidth', 3, 'DisplayName', 'Bathymetry');
ylabel('Bathymetry (m)', 'fontname', 'arial'); 
ylim([min(bathymetry(:,2))-1, 10]);

%%% 3.5: title and legend
% title('T = {T} s');
title(['T = ', num2str(T), ...
    ' s, h_{toe} = ', num2str(h_toe), ' m'],'fontname', 'arial');

ldg = legend([p1,p2,p3,p4,p5],'fontname', 'arial', 'Location', 'east');
ldg.FontSize = 18;

grid on
grid minor

%%% 3.6: Export figure
filename = 'fig1_wave_height_propagation.pdf'; % define the output file name

exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);

%% 4.0: Compute contributions from each wave dissipation mechanism

%%% 4.1: Compute each component
K_comp_veg = drag_components_f(aw_list_veg);
K_comp_noveg = drag_components_f(aw_list_noveg);

%%% 4.2: Plot mechanism breakdown for no vegetation case 
figure 

% vegetation types
p1 = area([bounds(1,1)-dx bounds(1,2)-dx], [-110 0], 'FaceColor', 'g', ...
    'FaceAlpha',0.1, 'EdgeColor', 'g', 'EdgeAlpha', 0.1,'DisplayName', veg_types{1});
hold on
area([bounds(1,1)-dx bounds(1,2)-dx], [0 110], 'FaceColor', 'g', ...
    'FaceAlpha',0.1, 'EdgeColor', 'g', 'EdgeAlpha', 0.1,'DisplayName', veg_types{1});
p2 = area([bounds(2,1)-dx bounds(2,2)-dx], [-110 0], 'FaceColor', 'g', ...
    'FaceAlpha',0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.3,'DisplayName', veg_types{2});
area([bounds(2,1)-dx bounds(2,2)-dx], [0 110], 'FaceColor', 'g', ...
    'FaceAlpha',0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.3,'DisplayName', veg_types{2});

% wave height
p3 = plot(aw_list_noveg(:,1), K_comp_noveg(:,2)*2, LineWidth=4, color = '#D95319', DisplayName='Shoaling');
p4 = plot(aw_list_noveg(:,1), K_comp_noveg(:,3)*2, LineWidth=4, color = '#EDB120', DisplayName='Wave breaking');
p5 = plot(aw_list_noveg(:,1), K_comp_noveg(:,4)*2, LineWidth=4, color = '#7E2F8E', DisplayName='Bed friction');
ylabel({'Contributed wave height change (m)'});
ylim([-0.4 0.2]);

title({'Mechanism breakdown for without vegetation'})
xlabel('Distance from seawall (m)', 'fontname', 'arial'); % wave height
set(gca,'FontSize', 20, 'color', 'none');

grid on;
grid minor;

legend([p1, p2, p3, p4, p5],{veg_types{1},veg_types{2},'Shoaling', 'Wave breaking','Bed friction'}, Location='southeast');

%%% 4.3: Export figure
filename = 'fig2_breakdown_without_veg.pdf'; % define the output file name

exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);


%%% 4.4: Plot mechanism breakdown for with vegetation case 
figure 

% vegetation types
p1 = area([bounds(1,1)-dx bounds(1,2)-dx], [-110 0], 'FaceColor', 'g', ...
    'FaceAlpha',0.1, 'EdgeColor', 'g', 'EdgeAlpha', 0.1,'DisplayName', veg_types{1});
hold on
area([bounds(1,1)-dx bounds(1,2)-dx], [0 110], 'FaceColor', 'g', ...
    'FaceAlpha',0.1, 'EdgeColor', 'g', 'EdgeAlpha', 0.1,'DisplayName', veg_types{1});
p2 = area([bounds(2,1)-dx bounds(2,2)-dx], [-110 0], 'FaceColor', 'g', ...
    'FaceAlpha',0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.3,'DisplayName', veg_types{2});
area([bounds(2,1)-dx bounds(2,2)-dx], [0 110], 'FaceColor', 'g', ...
    'FaceAlpha',0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.3,'DisplayName', veg_types{2});

% wave height:
p3 = plot(aw_list_veg(:,1), K_comp_veg(:,1)*2, LineWidth=4, color='#0072BD', DisplayName='Vegetation');
hold on
p4 = plot(aw_list_veg(:,1), K_comp_veg(:,2)*2, LineWidth=4, color = '#D95319', DisplayName='Shoaling');
p5 = plot(aw_list_veg(:,1), K_comp_veg(:,3)*2, LineWidth=4, color = '#EDB120', DisplayName='Wave breaking');
p6 = plot(aw_list_veg(:,1), K_comp_veg(:,4)*2, LineWidth=4, color = '#7E2F8E', DisplayName='Bed friction');
ylabel({'Contributed wave height change (m)'});
ylim([-1.0 0.2]);

title({'Mechanism breakdown for with vegetation'})
xlabel('Distance from seawall (m)', 'fontname', 'arial'); % wave height
set(gca,'FontSize', 20, 'color', 'none');

grid on;
grid minor;

legend([p1, p2, p3, p4, p5, p6],{veg_types{1},veg_types{2},'Vegetation','Shoaling', ...
    'Wave breaking','Bed friction'}, Location='southeast');

%%% 4.5: Export figure
filename = 'fig3_breakdown_with_veg.pdf'; % define the output file name
set(gcf, 'PaperPositionMode', 'auto');  % Ensure the paper position mode is auto
exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);