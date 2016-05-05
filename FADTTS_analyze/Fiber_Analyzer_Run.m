function [] = Fiber_Analyzer_Run(Fiber_Analyzer_Results_Directory, wNIC_groups_file)

% Fiber_Analyzer_Results_Directory:
% The "FiberProfile" directory that
% results from a run of DTI Fiber Atlas Analyzer

% wNIC_groups_file: the csv groups file for each case.
% This file would have one column that contains IDs for each case,
% corresponding to the first line of a DTI Fiber Profile; the first
% row of this CSV would be "id,<group_1>,<group_2>,..." and so on.
% The rest are binary indications of whether each case belongs to a
% certain group.

% Add the path to FADTTS, if it isn't already present
%addpath '/SAN/vision/camino5/camino/MattLeming/lib/FADTTS/';

% To make sure that functions aren't lost upon a change in directory.
addpath(pwd);

% The colors of different groups, to make everything numerical
color_map = {'r', 'b', 'g', 'y', 'm', 'k'};

% Measurements (currently this script is only tested with one measurement, 'FA', so that will need further testing
diffusion_properties = {'FA','RD','AD','MD'};

figure_number=1; % This number increases for every new graph
%threshold = 0.05; % For putting dots on certain graphs if they go below a p-value threshold

% Use 100 when testing scripts. Use 10000 when running scripts for real.
number_of_permutations = 100;

% Diffusion Properties
% number of different properties (parameters) to be analyzed (ex. FA,RD,AD,MD)
number_of_diffusion_properties = 1;
% if we consider m response variables, cells(m,1).
measurements=cell(number_of_diffusion_properties,1);

for i=1:number_of_diffusion_properties
    measurements{i}=diffusion_properties{i};
end

working_directory = dir(Fiber_Analyzer_Results_Directory);
working_directory = working_directory(arrayfun(@(x) x.name(1), working_directory) ~= '.');
cd(Fiber_Analyzer_Results_Directory)
for i = 1: length(working_directory)
    %cd (fullfile(Fiber_Analyzer_Results_Directory, working_directory(i).name));

    % name_select takes in a regular expression to help it identify the
    % fiber profile and the groups file in a certain directory, using that
    % to assign the fiber name as well. In this case, the group file is 
    % already assigned, but it can identify the fiber profile by
    % identifying "fa" at the beginning of the filename. See name_select
    % for specifics.
    [fiber_name,wNIC_file,~] = name_select('^fa', '');
    
    % Plots every case, color coded by group
    binary_group_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map);
    figure_number = figure_number + 1;

    %Untested
    %beta_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map)
    %figure_number = figure_number + 1;
    
    % Plots averages and standard deviations
    group_avg_std_dev_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map);
    figure_number = figure_number + 1;
    
    %global_p_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map, number_of_permutations);
    %figure_number = figure_number + 1;
    
    % Plots the corrected local p-values.
    corrected_p_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map, number_of_permutations);
    figure_number = figure_number + 1;
    
end
