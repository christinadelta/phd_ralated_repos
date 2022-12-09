function [pval slope, h] = myScatter(x,y,robust,mcol,mshape,linestyle)
% [pval slope] = myScatter(x,y,robust,mcol,mshape)
% plot a scatterplot, and regression line, and return p-value of the
% beta estimate of the slope
% rtfitweight: type of robust regression. Usually using cauchy. Options:
% 'andrews','bisquare','cauchy','fair','huber','logistic','talwar','welsch';
%
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2012-2014 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------

robusttype = 'cauchy';

if ~exist('robust'); robust = false; end
if ~exist('mcol'); mcol = [0 0 0]; end
if ~exist('mshape'); mshape = 'd'; end
if ~exist('linestyle'); linestyle = '-'; end

scatter(x,y,[],mcol,mshape);
hold on
if robust
    [B, stats] = robustfit(x,y,robusttype);
else
    [B,dev,stats] = glmfit(x,y);
end
pval = stats.p(2);
slope = B(2);

lx = [min(x) max(x)];
ly = lx*slope+B(1);
h = plot(lx,ly,linestyle,'linewidth',2,'color',mcol);

xlim([min(x)-range(x)/20 max(x)+range(x)/20]);
% 
% xl = xlim;
% yl = ylim;
% 
% tmp = [min(xl(1),yl(1)),max(xl(2),yl(2))];
% newlim = [(tmp(1) - diff(tmp)*0.1),(tmp(2) + diff(tmp)*0.1)];

txt = sprintf('p = %0.02f, slope = %0.02f',pval,slope);
% title(txt);
legend(txt,'Location','SouthOutside');
legend boxoff
% text(max(x) + range(x)/10,max(y) + range(y)/10,txt);

end
