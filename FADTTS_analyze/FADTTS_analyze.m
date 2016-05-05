function [] = FADTTS_analyze(target_directory, binary_groups, autoselect)

%addpath '/SAN/vision/camino5/camino/MattLeming/lib/FADTTS/';
addpath(pwd);
cd (target_directory);

color_map = {'-r', '-b', '-g', '-y', '-m', '-k'};
diffusion_properties = {'FA','RD','AD','MD'};

fignum=1;
threshold = 0.05;

% number of bootstrapping permutations to use in Hypothesis testing
numPerms = 500; % Use 100 when testing scripts. Use 10000 when running scripts for real.
myft=1; % defining my fiber tract variable to use in file names

ft=cell(1,1);


% Diffusion Properties
nofDiffProps = 1; % number of different properties (parameters) to be analyzed (ex. FA,RD,AD,MD)
Mnames=cell(nofDiffProps,1);  % if we consider m response variables, cells(m,1).
Mnames{1}='FA';
%Mnames{2}='RD';
%Mnames{3}='AD';
%Mnames{4}='MD';

if(autoselect)
  [ft{1},wNIC_file,wNIC_groups_file] = name_select('FA_wNIC', 'groups_wNIC');
end

if binary_groups
  binary_group_plot(wNIC_file, wNIC_groups_file, ft{1}, fignum, Mnames, color_map);
  fignum = fignum + 1;
end

beta_plot(wNIC_file, wNIC_groups_file, ft{1}, fignum, Mnames, color_map)
fignum = fignum + 1;

group_avg_std_dev_plot(wNIC_file, wNIC_groups_file, ft{1}, fignum, Mnames, color_map);
fignum = fignum + 1;

global_p_plot(wNIC_file, wNIC_groups_file, ft{1}, fignum, Mnames, color_map, numPerms);
fignum = fignum + 1;

corrected_p_plot(wNIC_file, wNIC_groups_file, ft{1}, fignum, Mnames, color_map, numPerms);
fignum = fignum + 1;




%%%%%%%%%%%

[params, myparams, mycovars, VarName, dataFiber1All, diffusionFiles, designdata] = read_fiber_data(wNIC_file,wNIC_groups_file, Mnames);
arclength = dataFiber1All(:,1); % take first column => arclength from dtitractstatCLP fiber file
tractdata=[arclength zeros(size(arclength,1),1) zeros(size(arclength,1),1)];
nofeatures=size(diffusionFiles, 1);
[NoSetup, arclength_allPos, Xdesign, Ydesign] = MVCM_read(tractdata, designdata, diffusionFiles, nofeatures);
NoSetupT = num2cell(NoSetup);
[n,L0,p,m] = deal(NoSetupT{:})
[mh]=MVCM_lpks_wob(NoSetup, arclength_allPos, Xdesign, Ydesign);
[efitBetas, efitBetas1, InvSigmats, efitYdesign]=MVCM_lpks_wb1(NoSetup, arclength_allPos, Xdesign, Ydesign, mh);
ResYdesign=Ydesign-efitYdesign;
[ResEtas,efitEtas,eSigEta]=MVCM_sif(arclength_allPos,ResYdesign); %Matt: replace 'arclength_allPos' with 'arclength'?
[mSigEtaEig, mSigEta]=MVCM_eigen(efitEtas);
[ebiasBetas] = MVCM_bias(NoSetup,arclength_allPos,Xdesign,Ydesign,InvSigmats,mh);
 for pp=2:p
    %individual and global statistics calculation
    cdesign=zeros(1,p);
    cdesign(pp)=1;
    Cdesign=kron(eye(m),cdesign);
    B0vector=zeros(m,L0);
    [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arclength_allPos,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
    Gstats(1,pp-1)=Gstat;
    Lstats(:,pp-1)=Lstat;
    % Generate random samples and calculate the corresponding statistics and pvalues
    numPerms;
    [Gpval] = MVCM_bstrp_pvalue3(NoSetup,arclength_allPos,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,numPerms); 
    Gpvals(1,pp-1)=Gpval;
 end
Lpvals=1-chi2cdf(Lstats,m);
Lpvals_FDR=zeros(size(Lpvals));
for i=1:mycovars
  Lpvals_FDR(:,i)=mafdr(Lpvals(:,i),'BHFDR',true);
end
%% Omnibus Covariate Confidence Bands

[Gvalue] = MVCM_cb_Gval(arclength_allPos,Xdesign,ResYdesign,InvSigmats,mh,numPerms);
alpha=0.05;
[CBands] = MVCM_CBands(n,alpha,Gvalue,efitBetas,zeros(size(ebiasBetas))); % new Conf Bands formula

%% plot omnibus confidence bands

ylabs=cell(1+mycovars,1);
ylabs{1}=sprintf('95 percent confidence band for %s intercept',ft{myft});     

%
%% 4. Post-hoc Hypothesis Test --> Which Diffusion Parameter is significant where for each covariate? A univariate test in a multivariate model.

% for comparing the significance of each diffusion parameter for each covariate
posthoc_Gpvals=zeros(m,p-1);
posthoc_Lpvals=zeros(L0,m,p-1);

for pii = 2:p; % each covariate's betas
    for mii=1:m
        Cdesign=zeros(1,m*p);
        Cdesign(1+(mii-1)*p+(pii-1))=1;
        B0vector=zeros(1,L0);
        [Gstat,Lstat] = MVCM_ht_stat(NoSetup,arclength_allPos,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
        % Generate random samples and calculate the corresponding statistics and pvalues
        GG=numPerms;
        posthoc_Gpvals(mii,pii-1)= MVCM_bstrp_pvalue3(NoSetup,arclength_allPos,Xdesign,Ydesign,efitBetas1,InvSigmats,mh,Cdesign,B0vector,Gstat,GG);
        posthoc_Lpvals(:,mii,pii-1)=1-chi2cdf(Lstat,1);
    end
end

% Global P values for posthoc test

posthoc_Gpvals % for FA, RD for each covariate 

% Save Post-hoc test Global p-values for each diffusion parameter

csvwrite(sprintf('%s_%s_%s_posthoc_Global_pvalues.csv',ft{myft},params{myparams},VarName{1}),posthoc_Gpvals); %save csv file 

%% correct posthoc test local p-values with FDR 
% this corrects the posthoc local p-values for multiple comparisons

posthoc_Lpvals_FDR=zeros(size(posthoc_Lpvals));
for mii=1:m
    for pii=1:mycovars
        posthoc_Lpvals_FDR(:,mii,pii)=mafdr(posthoc_Lpvals(:,mii,pii),'BHFDR',true);
    end
end

% save Corrected Local P-Values csv file
for mii=1:m
    csvwrite(sprintf('%s_%s_posthoc_Local_pvalues.csv',ft{myft},Mnames{mii}),posthoc_Lpvals(:,mii,:)); 
    csvwrite(sprintf('%s_%s_posthoc_corrected_Local_pvalues.csv',ft{myft},Mnames{mii}),posthoc_Lpvals_FDR(:,mii,:)); 
end

%% Plot CORRECTED Post-hoc test Local P-values for Each Difffusion Parameter

for pii = 2:p; % each covariate's betas
    figure(fignum);
    for mii=1:m
    plot(arclength,-log10(posthoc_Lpvals_FDR(:,1,pii-1)),strcat(color_map{mii}, '.'),'LineWidth', 2,'MarkerSize',15); % FA
    end
    hold off;
    xL = get(gca,'XLim');
    line(xL,[1.3 1.3],'Color','black'); % line at 1.3 to mark significance level
    quickplot(h, posthoc_Lpvals_FDR, 'arclength', '-log10(p)', arclength, [Mnames(1:nofDiffProps)], sprintf('%s_%s_%s_posthoc_corrected_Local_pvalues.pdf',ft{myft},params{myparams},char(VarName(pii-1))));
end

 
 %% plot CORRECTED Post-Hoc Test SIGNIFICANT BETAS only- ZOOOMED IN - for ALL diffusion paramters on one plot -- A plot for each COVARIATE

