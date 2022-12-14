function simout = analysis_wsls(a, r)

% the 
alast       = [nan a(1:end-1)];
stay        = alast == a;
rlast       = [nan r(1:end-1)];

winstay     = nanmean(stay(rlast==1));
losestay    = nanmean(stay(rlast==0));

simout      = [losestay winstay];




end % end of fuction 