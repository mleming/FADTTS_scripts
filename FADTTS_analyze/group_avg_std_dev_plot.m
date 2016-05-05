function [] = group_avg_std_dev_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map)
if ~exist('wNIC_file', 'var')
    [~,wNIC_file,~] = name_select('FA', 'groups');
end

if ~exist('wNIC_groups_file', 'var')
    [~,~,wNIC_groups_file] = name_select('FA', 'groups');
end

%if ~exist('fiber_name', 'var')
%    [fiber_name,~,~] = name_select('FA', 'groups');
%end

if ~exist('measurements', 'var')
    measurements=cell(1,1);  % if we consider m response variables, cells(m,1).
    measurements{1}='FA';
end

if ~exist('figure_number', 'var')
    figure_number = 1;
end

if ~exist('color_map', 'var')
    color_map = {'r', 'b', 'g', 'y', 'm', 'k'};
end

[~, ~, ~, VarName, all_data_fibers, diffusion_files, design_data] = read_fiber_data(wNIC_file,wNIC_groups_file,measurements);
arc_index_points = all_data_fibers(:,1); % take first column => arclength from dtitractstatCLP fiber file
tract_data=[arc_index_points zeros(size(arc_index_points,1),1) zeros(size(arc_index_points,1),1)];
[~, ~, ~, Ydesign] = MVCM_read(tract_data, design_data, diffusion_files, size(diffusion_files, 1));

all_average = [];
all_standard_deviation = [];


figure(figure_number);


for mii=1:size(measurements,1)
    for i = 2: size(design_data,2)
        all_average = [all_average mean(Ydesign(design_data(:,i)==1,:,mii))'];
        all_standard_deviation = [all_standard_deviation std(Ydesign(design_data(:,i)==1,:,mii))'];
    end
    for i = 1:length(all_average(mii, 1:end))
        temp_average = all_average(1:end, i);
        temp_standard_deviation = all_standard_deviation(1:end, i);
        h(i)=plot(arc_index_points, temp_average(:,:,mii),strcat('-',color_map{i}),'LineWidth', 2); % Cocaine average
        hold on;
    end
    for i = 1:length(all_average(mii, 1:end))
        temp_average = all_average(1:end, i);
        temp_standard_deviation = all_standard_deviation(1:end, i);
        plot(arc_index_points, temp_average(:,:,mii)+temp_standard_deviation(:,:,mii),strcat('--',color_map{i}),'LineWidth', 2);
        hold on;
        plot(arc_index_points, temp_average(:,:,mii)-temp_standard_deviation(:,:,mii),strcat('--',color_map{i}),'LineWidth', 2);
        hold on;
    end
    hold off;
    quickplot(h, [arc_index_points all_average], 'arclength', measurements{mii}, arc_index_points, VarName, sprintf('%s %s Group Average and Standard Deviation 2group',fiber_name,measurements{mii}));
    figure_number = figure_number + 1;
end