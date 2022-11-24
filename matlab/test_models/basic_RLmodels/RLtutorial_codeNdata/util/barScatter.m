function [hb,he]=barScatter(data, errors, bw_xtick, bScatter)
% [hb,he]=barScatter(data, errors, bw_xtick, bScatter)
%
% barScatter is the m-by-n matrix of data to be plotted. barScatter calls
% the MATLAB bar function and plots
% - m bars
% - a scatterplot of the individual datapoints, if bScatter = true
% - errorbars of length error (vector of m long), or, if empty, using SD
% - xtick labels on the x axis. 
% NB note that the grouping occurs by row, so each row is a separate group
%
% See the MATLAB functions bar and errorbar for more information
%
% [hb,he]=barScatter(...) will give handle HB to the bars and HE to the
% errorbars.
% 
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------

% Get function arguments
if nargin < 1
    error('Must have at least the first argument:  barweb(data, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend)');
elseif nargin == 1
    bw_xtick = [];
    bScatter = true;
elseif nargin == 2
    bw_xtick = [];
    bScatter = true;
elseif nargin == 3
    bScatter = true;
end

barvalues = squeeze(nanmean(data,2));
nbars = size(data,1);

% Plot bars and errors
hb=bar(barvalues);
hold on;

set(gca, 'fontsize', 10, 'fontweight', 'bold');
he=[];
if isempty(errors) | nargin==1;
    for t = 1:nbars
        tmp = data(t,~isnan(data(t,:)));
        errors(t) = std(tmp);
%         errors(t) = errors(t)./sqrt(length(tmp)); % compute standard error of the mean
    end
elseif length(barvalues) ~= length(errors)
    error('barvalues and errors vectors must be of same dimension');
end

if bScatter
    for t = 1:nbars
        tmp = data(t,~isnan(data(t,:)));
        s = scatter(t*ones(1,length(tmp)),tmp,[],'k', 'jitter','on', 'jitterAmount',0.05);
        set(s,'MarkerEdgeColor',[1 0 0],'linewidth',2);
    end
end
th=errorbar(barvalues, errors, 'k', 'linestyle', 'none','linewidth',2.5);
set(th,'Linewidth',2.5,'color',[.4 .4 .4]);
he=[he;th];
removeErrorBarEnds(th)
xlim([0.5 nbars+0.5]);

if~isempty(bw_xtick)
    set(gca,'xticklabel', bw_xtick);
end

hold off;

return
