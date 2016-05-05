function [] = global_p_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map, numPerms)
  if ~exist('fiber_name', 'var')
    [fiber_name,o,f] = name_select('FA_wNIC', 'groups_wNIC');
  end

  if ~exist('measurements', 'var')
    measurements=cell(1,1);  % if we consider m response variables, cells(m,1).
    measurements{1}='FA';
  end
    
  if ~exist('figure_number', 'var')
    figure_number = 1;
  end
    
  if ~exist('color_map', 'var')
    color_map = {'r', 'b', 'g', 'y', 'm', 'k','-r', '-b', '-g', '-y', '-m', '-k'};
  end
  
  if ~exist('numPerms', 'var')
    numPerms = 500;
  end
  
  [params, myparams, ~, variable_names, all_fiber_data, diffusion_files, designdata] = read_fiber_data(wNIC_file,wNIC_groups_file, measurements);
  arc_index_points = all_fiber_data(:,1); % take first column => arclength from dtitractstatCLP fiber file
  tract_data=[arc_index_points zeros(size(arc_index_points,1),1) zeros(size(arc_index_points,1),1)];
  [NoSetup, arc_index_points_all_positions, Xdesign, Ydesign] = MVCM_read(tract_data, designdata, diffusion_files, size(diffusion_files, 1));
  NoSetupT = num2cell(NoSetup);
  [~,L0,p,m] = deal(NoSetupT{:});
  [mh]=MVCM_lpks_wob(NoSetup, arc_index_points_all_positions, Xdesign, Ydesign);
  [efitBetas, efitBetas1, InvSigmats, efitYdesign]=MVCM_lpks_wb1(NoSetup, arc_index_points_all_positions, Xdesign, Ydesign, mh);
  [~,~,eSigEta]=MVCM_sif(arc_index_points_all_positions,Ydesign-efitYdesign);
  [ebiasBetas] = MVCM_bias(NoSetup,arc_index_points_all_positions,Xdesign,Ydesign,InvSigmats,mh);

  Lstats=zeros(L0,p-1);
  Gpvals=zeros(1,p-1);
  
  for i=2:p
    %individual and global statistics calculation
    cdesign=zeros(1,p);
    cdesign(i)=1;
    Cdesign=kron(eye(m),cdesign);
    B0vector=zeros(m,L0);
    [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arc_index_points_all_positions,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
    Lstats(:,i-1)=Lstat;

    % Generate random samples and calculate the corresponding statistics and pvalues
    [Gpval] = MVCM_bstrp_pvalue3(NoSetup,arc_index_points_all_positions,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,numPerms); 
    Gpvals(1,i-1)=Gpval;
  end
  Lpvals=1-chi2cdf(Lstats,m);
  figure(figure_number);
  xL = get(gca,'XLim');
  line(xL,[1.3 1.3],'Color','black');  % line at 1.3 to mark significance level
  hold on;
  for i = 1:length(Lpvals(1,1:end))
        h(i) = plot(arc_index_points,-log10(Lpvals(:,i)),strcat('-',color_map{i},'.'),'LineWidth', 2,'MarkerSize',15); % Group Cocaine
        hold on;
  end 
  %h(1) = plot(arc_index_points,-log10(Lpvals(:,1)),'-b.','LineWidth', 2,'MarkerSize',15); % Group Cocaine
  hold off;
  quickplot(h, Lpvals, 'arclength', '-log10(p)', arc_index_points, variable_names', sprintf('%s %s Local p-values',fiber_name,params{myparams}));
  csvwrite(sprintf('%s_%s_%s_Global_pvalues.csv',fiber_name,params{myparams},variable_names{1}),Gpvals);
 % figure_number = figure_number + 1;
 % figure(figure_number);
 % quickplot(h, Gpvals, 'arclength', '-log10(p)', arc_index_points, variable_names', sprintf('%s %s Global p-values',fiber_name,params{myparams}));      
     
