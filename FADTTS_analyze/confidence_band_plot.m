function o = confidence_band_plot(wNIC_file, wNIC_groups_file, fiber_name, figure_number, measurements, color_map, number_of_permutations)
if ~exist('fiber_name', 'var')
    [fiber_name,~,~] = name_select('FA_wNIC', 'groups_wNIC');
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

if ~exist('number_of_permutations', 'var')
    number_of_permutations = 500;
end

[~, ~, number_of_covariates, variable_names, dataFiber1All, diffusionFiles, designdata] = read_fiber_data(wNIC_file,wNIC_groups_file, measurements);
arclength = dataFiber1All(:,1); % take first column => arclength from dtitractstatCLP fiber file
tractdata=[arclength zeros(size(arclength,1),1) zeros(size(arclength,1),1)];
[NoSetup, arclength_allPos, Xdesign, Ydesign] = MVCM_read(tractdata, designdata, diffusionFiles, size(diffusionFiles, 1));
NoSetupT = num2cell(NoSetup);
[n,L0,p,~] = deal(NoSetupT{:});
[mh]=MVCM_lpks_wob(NoSetup, arclength_allPos, Xdesign, Ydesign);
[efitBetas, ~, InvSigmats, efitYdesign]=MVCM_lpks_wb1(NoSetup, arclength_allPos, Xdesign, Ydesign, mh);
[~,~,eSigEta]=MVCM_sif(arclength_allPos,Ydesign-efitYdesign); %Matt: replace 'arclength_allPos' with 'arclength'?
[ebiasBetas] = MVCM_bias(NoSetup,arclength_allPos,Xdesign,Ydesign,InvSigmats,mh);

Lstats=zeros(L0,p-1);
for pi=2:p
    %individual and global statistics calculation
    cdesign=zeros(1,p);
    cdesign(pi)=1;
    Cdesign=kron(eye(size(measurements,1)),cdesign);
    B0vector=zeros(size(measurements,1),L0);
    [~,Lstat] = MVCM_ht_stat(NoSetup,arclength_allPos,Xdesign,efitBetas,eSigEta,Cdesign,B0vector,ebiasBetas);
    Lstats(:,pi-1)=Lstat;
end
local_p_values=1-chi2cdf(Lstats,size(measurements,1));
local_p_values_FDR=zeros(size(local_p_values));
for i=1:number_of_covariates
    local_p_values_FDR(:,i)=mafdr(local_p_values(:,i),'BHFDR',true);
end

%% Omnibus Covariate Confidence Bands

[Gvalue] = MVCM_cb_Gval(arclength_allPos,Xdesign,Ydesign-efitYdesign,InvSigmats,mh,number_of_permutations);
alpha=0.05;
[CBands] = MVCM_CBands(n,alpha,Gvalue,efitBetas,zeros(size(ebiasBetas))); % new Conf Bands formula

%% plot omnibus confidence bands

y_axis_labels=cell(1+number_of_covariates,1);
y_axis_labels{1}=sprintf('95 percent confidence band for %s intercept',fiber_name);
for i=2:(number_of_covariates+1)
    y_axis_labels{i}=sprintf('95 percent confidence band for %s %s', fiber_name, char(variable_names(i-1)));
end
for mi=1:size(measurements,1)
    for pi=1:p
        figure(figure_number);
        figure_number = figure_number + 1;
        plot(arclength,efitBetas(pi,:,mi),strcat('-',color_map{2}),'LineWidth', 2); % Parameter (FA,RD)
        hold on;
        plot(arclength,CBands(2*pi-1,:,mi),strcat('--',color_map{1}),'LineWidth', 2); % Lower Conf Band
        hold on;
        plot(arclength,CBands(2*pi,:,mi),strcat('--',color_map{1}),'LineWidth', 2); % Upper Conf Band
        xlabel('arclength');
        ylabel(measurements{mi});
        title(y_axis_labels{pi},'Interpreter','none');
        xlim([min(arclength) max(arclength)]);
        xL = get(gca,'XLim');
        line(xL,[0 0],'Color','black'); % line at zero
        hold off;
        % save plot
        if (pi==1)
            figurename=sprintf('%s_%s_%s_confidence_band.pdf',fiber_name,measurements{mi},'intercept');
        else
            figurename=sprintf('%s_%s_%s_confidence_band.pdf',fiber_name,measurements{mi}, char(variable_names(pi-1)));
        end
        saveas(gcf,figurename,'pdf');
    end
end
