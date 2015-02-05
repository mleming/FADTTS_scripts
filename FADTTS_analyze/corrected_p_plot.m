function o = corrected_p_plot(wNIC_file, wNIC_groups_file, fiber_name, fignum, measurements, color_map, number_of_permutations)


if ~exist('wNIC_file', 'var')
    [~,wNIC_file,~] = name_select('fa', 'groups');
end

if ~exist('wNIC_groups_file', 'var')
    [~,~,wNIC_groups_file] = name_select('fa', 'groups');
end

if ~exist('fiber_name', 'var')
    [fiber_name,~,~] = name_select('fa', 'groups');
end

if ~exist('measurements', 'var')
    measurements=cell(1,1);  % if we consider m response variables, cells(m,1).
    measurements{1}='FA';
end

if ~exist('fignum', 'var')
    fignum = 1;
end

if ~exist('color_map', 'var')
    color_map = {'r', 'b', 'g', 'y', 'm', 'k'};
end

if ~exist('number_of_permutations', 'var')
    number_of_permutations = 500;
end




[params, myparams, my_covariates, variable_names, all_fiber_data, diffusion_files, design_data] = read_fiber_data(wNIC_file,wNIC_groups_file, measurements);
arclength = all_fiber_data(:,1); % take first column => arclength from dtitractstatCLP fiber file
tractdata=[arclength zeros(size(arclength,1),1) zeros(size(arclength,1),1)];
nofeatures=size(diffusion_files, 1);
[NoSetup, arclength_allPos, Xdesign, Ydesign] = MVCM_read(tractdata, design_data, diffusion_files, nofeatures);
NoSetupT = num2cell(NoSetup);
[n,L0,p,m] = deal(NoSetupT{:});
[mh]=MVCM_lpks_wob(NoSetup, arclength_allPos, Xdesign, Ydesign);
[efitBetas, efitBetas1, InvSigmats, efitYdesign]=MVCM_lpks_wb1(NoSetup, arclength_allPos, Xdesign, Ydesign, mh);
ResYdesign=Ydesign-efitYdesign;
[ResEtas,efitEtas,eSigEta]=MVCM_sif(arclength_allPos,ResYdesign); %Matt: replace 'arclength_allPos' with 'arclength'?
[ebiasBetas] = MVCM_bias(NoSetup,arclength_allPos,Xdesign,Ydesign,InvSigmats,mh);
[mSigEtaEig, mSigEta]=MVCM_eigen(efitEtas);

Gstats=zeros(1,p-1);
Lstats=zeros(L0,p-1);
Gpvals=zeros(1,p-1);

for pi=2:p
    %individual and global statistics calculation
    cdesign=zeros(1,p);
    cdesign(pi)=1;
    Cdesign=kron(eye(m),cdesign);
    B0vector=zeros(m,L0);
    [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arclength_allPos,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
    Gstats(1,pi-1)=Gstat;
    Lstats(:,pi-1)=Lstat;
    %[Gpval] = MVCM_bstrp_pvalue3(NoSetup,arclength_allPos,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,number_of_permutations);   
    %Gpvals(1,pi-1)=Gpval;
    % Generate random samples and calculate the corresponding statistics and pvalues
end
local_p_values=1-chi2cdf(Lstats,m);
local_p_values_FDR=zeros(size(local_p_values));
for i=1:my_covariates
    local_p_values_FDR(:,i)=mafdr(local_p_values(:,i),'BHFDR',true);
end
global_p_values_FDR=zeros(size(Gpvals));
for i=1:size(Gpvals)
   global_p_values_FDR(:,i)=mafdr(Gpvals(:,i), 'BHFDR', true);
end

%% plot corrected local p values for each covariate
% need a plot line for each covariate
figure(fignum);
xL = get(gca,'XLim');
line(xL,[1.3 1.3],'Color','black'); % line at 1.3 to mark significance level
line(xL,[1.3 1.3],'Color','black');
hold on;
    for pi=2:p
    h(pi-1) = plot(arclength,-log10(local_p_values_FDR(:,pi-1)),strcat('-',color_map{pi-1},'.'),'LineWidth', 2,'MarkerSize',15);
    hold on;
    end
hold off;
quickplot(h, local_p_values_FDR, 'arclength', '-log10(p)', arclength, variable_names', sprintf('%s %s Corrected Local p-values 2group',fiber_name,params{myparams}));

%fignum = fignum + 1;
%figure(fignum);
%line(xL,[1.3 1.3],'Color','black');
%hold on;
%h(1) = plot(arclength, -log10(global_p_values_FDR'),  strcat('-',color_map{1},'.'), 'LineWidth', 2, 'MarkerSize', 15);
%hold off;
%quickplot(h, local_p_values_FDR, 'arclength', '-log10(p)', arclength, cat(1, 'Intercept', variable_names'), sprintf('%s %s Corrected Global p-values',fiber_name,params{myparams}));








