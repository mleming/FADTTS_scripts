function [] = estimated_coefficient_plot(wNIC_file, wNIC_groups_file, fiber_name, fignum, Mnames, color_map, numPerms)
  if ~exist('fiber_name', 'var')
    [fiber_name,o,f] = name_select('FA_wNIC', 'groups_wNIC');
  end

  if ~exist('Mnames', 'var')
    Mnames=cell(1,1);  % if we consider m response variables, cells(m,1).
    Mnames{1}='FA';
  end
    
  if ~exist('fignum', 'var')
    fignum = 1;
  end
    
  if ~exist('color_map', 'var')
    color_map = {'-r', '-b', '-g', '-y', '-m', '-k'};
  end
  
  if ~exist('numPerms', 'var')
    numPerms = 500;
  end
  
  
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


%% plot ALL betas zoomed in - for a multivariate analysis - NEW way of plotting betas 
% betas are the coefficients that describe how related the covariate is to
% the parameter


 for mii=1:m
     beta_plot_new(arclength, efitBetas, Lpvals_FDR, threshold, mii, p, fignum, color_map);
     fignum = fignum + 1;   
     quickplot(h, efitBetas, 'arclength', Mnames{mii}, arclength, VarName, sprintf('%s %s %s Estimated Coefficient',ft{myft},Mnames{mii},VarName{mycovars}));
 end
