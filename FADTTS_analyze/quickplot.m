function o = quickplot(h, data, xl, yl, x_minmax, vname, figuretitle)
    xlabel(xl);
    ylabel(yl);
    xlim([min(x_minmax) max(x_minmax)]);
    title(figuretitle, 'Interpreter','none');
    legend(h, vname,'Location','Best');  
    % save plot
    csvwrite(strcat(figuretitle, '.csv'), data)
    saveas(gcf,strcat(figuretitle,'.pdf'),'pdf');
    clear h;