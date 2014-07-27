function o = beta_plot(wNIC_file, wNIC_groups_file, fiber_name, fignum, Mnames, color_map)

    if ~exist('fiber_name', 'var')
      [fiber_name,~,~] = name_select('FA_wNIC', 'groups_wNIC');
    end

    if ~exist('Mnames', 'var')
      Mnames=cell(1,1);  % if we consider m response variables, cells(m,1).
      Mnames{1}='FA';
    end
    
    if ~exist('fignum', 'var')
      fignum = 1;
    end
    
    if ~exist('color_map', 'var')
      color_map = {'r', 'b', 'g', 'y', 'm', 'k'};
    end
    
    [~, ~, ~, variable_names, dataFiber1All, diffusionFiles, designdata] = read_fiber_data(wNIC_file,wNIC_groups_file,Mnames);
    arclength = dataFiber1All(:,1); % take first column => arclength from dtitractstatCLP fiber file
    tractdata=[arclength zeros(size(arclength,1),1) zeros(size(arclength,1),1)];
    nofeatures=size(diffusionFiles, 1);
    [NoSetup, arclength_allPos, Xdesign, Ydesign] = MVCM_read(tractdata, designdata, diffusionFiles, nofeatures);
    NoSetupT = num2cell(NoSetup);
    [~,~,~,m] = deal(NoSetupT{:});
    [mh]=MVCM_lpks_wob(NoSetup, arclength_allPos, Xdesign, Ydesign);
    [efitBetas, ~, ~, ~]=MVCM_lpks_wb1(NoSetup, arclength_allPos, Xdesign, Ydesign, mh);

    
    
for mii=1:m
    figure(fignum);
    fignum = fignum + 1;
    h(1) = plot(arclength, efitBetas(1,:,mii),'-m','LineWidth', 2);
    for i=2:size(efitBetas,1)
    hold on
    h(i) = plot(arclength, efitBetas(i,:,mii),'-b','LineWidth', 2);
    end
    quickplot(h, efitBetas(:,:,mii), 'arclength', Mnames{mii}, arclength, cat(1, 'Intercept', variable_names') ,sprintf('%s %s Estimated Coefficients',fiber_name,Mnames{mii}));
    hold off;
end
