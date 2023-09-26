function h = plotAL(mean_lrs)

% plot mean learning rates for each stc level and volailtity 
% ----------
stable  = mean_lrs(:,1)';
vol     = mean_lrs(:,2)';
lrs     = [stable; vol];
colours = [[171,5,32]/256; [7, 104, 115]/256;];

h       = bar(lrs);
set(h,'FaceColor',colours(1,:)) % change colour of bars to red

end % end of function

