function hfig = getNewFig(varargin)
% hfig = getNewFig(varargin)
% arg1 = hfig, returns figure handle hfig+1. If empty, return handle 1
% arg2 = cmap
% -------------------------------------------------------------------------


if numel(varargin)<1
    hfig = 1;
else
hfig = varargin{1}+1;
end
hfig = figure; hold on; box off;

% paperSize = [10 60 500 500];
paperSize = [10 60 700 400];

if numel(varargin)==2
    cmap = varargin{2};
    set(hfig,'DefaultAxesColorOrder',cmap,'position',paperSize,...
        'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
else
    set(hfig,'position',paperSize,...
        'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
end

end