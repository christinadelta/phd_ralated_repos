function [xx, dp, trialstable, trialvolatile] = responseModel_v1(xstate, val)

% Created May 2023
% The function generates responses for the action-learning simulated data 
% to be used for modelling 

% ---------------

ii50 = xstate==0.5; 
xstate(ii50) = [];
dv = [0; val(1:end-1)]-.5;
p = 1./(1+exp(-dv));
p(ii50) = [];

N = length(p);
n = 10;

change_points = find(diff(xstate)~=0);
trialvolatile = false(N,1);
for i=1:length(change_points)
    trialvolatile(change_points(i) + (1:n)) = 1;    
end
trialstable = ~trialvolatile;
trialstable(1:n) = 0;

corr_action = xstate>=.5;

choice = p>=.5;
perf = choice==corr_action;


mpvol = mean(perf(trialvolatile,:));
mpstab = mean(perf(trialstable,:));

mpvol = mean(mpvol);
mpstab = mean(mpstab);

dp = mpstab - mpvol;

xx = [mpstab mpvol];

end % end of function