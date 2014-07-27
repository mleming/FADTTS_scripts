function [h, data] = beta_plot_new(arclength, efitBetas, Lpvals_FDR, threshold, mii, p, fignum, color_map)

figure(fignum);

%line
for pii=2:p
  plot(arclength,efitBetas(pii,:,mii),color_map{pii}); 
  hold on;
end

%significance circle
for pii = 2:p
  ind=find(Lpvals_FDR(:,pii-1)<=threshold);
  plot(arclength(ind),efitBetas(pii,ind,mii),strcat('o', color_map{pii}(2)),'LineWidth',2)
  hold on; 
end

%significance fill
for pii=2:p
      ind=find(Lpvals_FDR(:,pii-1)<=threshold);
       plot(arclength(ind),efitBetas(pii,ind,mii),strcat('*', color_map{pii}(2)),'LineWidth',2)
       hold on; 
end

xL = get(gca,'XLim');
line(xL,[0 0],'Color','black'); %line at zero

hold off;