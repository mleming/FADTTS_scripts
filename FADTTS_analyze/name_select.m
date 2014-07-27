function [fiber_name,wNIC_file,wNIC_groups_file] = name_select(wNIC_suffix, wNIC_groups_suffix)

directory = dir;
for i=1:length(directory)
    if (not(isempty(regexp(directory(i).name,wNIC_groups_suffix,'ONCE'))))
        wNIC_groups_file=directory(i).name;
    end
    
    if (not(isempty(regexp(directory(i).name,wNIC_suffix,'ONCE'))))
        wNIC_file=directory(i).name;
    end
end


if(~exist('wNIC_file', 'var') && ~exist('wNIC_groups_file', 'var'))
    error('Exit - files not found - check suffixes');
elseif(~exist('wNIC_file', 'var'))
    fiber_name = regexprep(regexprep(wNIC_groups_file, wNIC_groups_suffix,''),'[_\- ]*(.*)[_\- ]*\..*$','$1');
    wNIC_file = '';
elseif(~exist('wNIC_groups_file', 'var'))
    fiber_name = regexprep(regexprep(wNIC_file, wNIC_suffix, ''),'[_\- ]*(.*)[_\- ]*\..*$','$1');
    wNIC_groups_file = '';
else
    shortlen = length(wNIC_file);
    if(length(wNIC_groups_file) < shortlen)
        shortlen=length(wNIC_groups_file);
    end
    notfg = 1;
    for i=1:shortlen
        if (not(wNIC_file(i) == wNIC_groups_file(i)) && notfg == 1)
            notfg = 0;
            if(wNIC_file(i-1) == '_')
                fiber_name = wNIC_file(1:i-2);
            else
                fiber_name = wNIC_file(1:i-1);
            end
        end
    end
end
