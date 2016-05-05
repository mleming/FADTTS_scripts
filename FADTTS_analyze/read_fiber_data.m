function [params, myparams, number_of_covariates, covariates, all_fiber_data, diffusion_files, designdata]=read_fiber_data(wNIC_FA_file,wNIC_groups_file, measurement_names);

% Diffusion Properties Tested for file names  
params=cell(1,1);
params{1}='FA';
myparams=1; 

%% 1. we need to define 'designdata', 'diffusionFiles', and 'arclength'
% deleted 'tractdata'because this variable redfined the arclength using xyz
% coordinates from FiberViewer that we do not have. Instead, we use our
% arclength from dtitractstatCLP and shift it to start at zero.

%% COVARIATES
% this file has all covariates organized in columns:
% 1. GroupCocaine (CTL = 0, PCE = 1)
% 2. GroupNicotine (CTL = 0, Nic = 1)

%Covariates Tested for file names
group_file_id = fopen(wNIC_groups_file);
first_line = fgetl(group_file_id);
fclose(group_file_id);

delimiter='';
if (length(find(first_line=='	')) < length(find(first_line==',')))
 delimiter = ',';
else
  delimiter = '	';
end

number_of_covariates = length(find(first_line==delimiter));

diffusion_files=cell(size(measurement_names, 1),1);	% if we consider m properties (parameters: FA,MD,RD,AD), cells(m,1)

groups = read_mixed_csv(wNIC_groups_file, delimiter);
data = read_mixed_csv(wNIC_FA_file, delimiter);
covariates = groups(1,2:end);
remove = [];
for i=2:size(groups,1)
    if sum(ismember(data,groups(i,1))) == 0
        remove = [remove i];
    end
end
groups = removerows(groups,'ind',remove);
groups(2:end,:) = sortrows(groups(2:end,:));
data(:,2:end) = sortrows(data(:,2:end)')';
diffusion_files{1} = cellfun(@str2double, data(2:end, find(ismember(data(1,1:end), groups(2:end,1)))));
all_fiber_data = cellfun(@str2double, [data(2:end,1) data(2:end, find(ismember(data(1,1:end), groups(2:end,1))))]);
designdata = [ones(size(diffusion_files{1}(1,:),2),1) cellfun(@str2double, groups(2:end, 2:end))];



