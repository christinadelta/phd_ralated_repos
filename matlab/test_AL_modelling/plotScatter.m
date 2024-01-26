function f  = plotScatter(y, x, ylm, colours, param_title)

% for each block plot simulated and estimated parameter values
% y = estimated (fitted) parameter
% x = simulated parameter 

%% define initial parameters
rows    = 2;
cols    = 3;
xlm     = ylm;
clf
y = y'; x = x';
%% loop over groups (that is over blocks)
groups = size(y,1);

for ii = 1:groups

    % plot the easy and diff conditions seperately
    ax(ii)          = subplot(rows,cols,ii);

    yg              = y(ii,:);  % get mean across reps/subjects and transpose so that the results is a 1xn array
    xg              = x;        % same simulated x for each block 

    % add scatter 
    f               = scatter(ax(ii),xg,yg,colours,'*');
    
    % add line of fit 
    lfit(ii)        = lsline(ax(ii));
    lfit(ii).Color  = colours;
    lfit(ii).LineWidth  = 0.7;

    % axis labels
    ylabel('estimated parameter values')
    xlabel('simulated parameter values')

    tmp=corrcoef(xg,yg);
    str=sprintf('r= %1.2f',tmp(1,2));
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    text(-0.03,1.12, sprintf('block %d',ii),'Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?
    fontsize(gcf,12,"points")

        
end % end of groups loop

% title handle
str = "Parameter Recovery for ";
str2 = ' full RBPF';
t = sgtitle(append(str, param_title, str2));
t.FontSize = 20;


end % end of function 