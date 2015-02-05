function o = binary_group_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map)


if ~exist('wNIC_file', 'var')
    [~,wNIC_file,~] = name_select('fa', 'groups');
end

if ~exist('wNIC_groups_file', 'var')
    [~,~,wNIC_groups_file] = name_select('FA', 'groups');
end

if ~exist('fiber_name', 'var')
    [fiber_name,~,~] = name_select('FA', 'groups');
end

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

[~, ~, number_of_covariates, covariates, all_fiber_data, diffusion_files, design_data] = read_fiber_data(wNIC_file,wNIC_groups_file,measurements);
arclength = all_fiber_data(:,1); % take first column => arclength from dtitractstatCLP fiber file
tract_data=[arclength zeros(size(arclength,1),1) zeros(size(arclength,1),1)];
number_of_features=size(diffusion_files, 1);
[NoSetup, ~, ~, Ydesign] = MVCM_read(tract_data, design_data, diffusion_files, number_of_features);
temp = num2cell(NoSetup);
[n,~,~,~] = deal(temp{:});

for mii=1:size(measurements,1)
    figure(figure_number);
    data = design_data(:,2:end);
    for nii=1:n %subjects
        for dii = 1:number_of_covariates
            if (data(nii, dii) > 0)
                h(dii) = plot(arclength,Ydesign(nii,:,mii),strcat('-',color_map{dii}),'LineWidth', 2);
                hold on;
            end
        end
    end
    quickplot(h, data, 'arclength', measurements{mii}, arclength, covariates(1:number_of_covariates),sprintf('%s %s 2group',fiber_name,measurements{mii}));
    clear h;
    figure_number = figure_number + 1;
end