function [xx, dp, trialstable, trialvolatile,p] = responseModel_v1(xstate, val,tvolatile, tstable)

% Created May 2023
% The function generates responses for the action-learning simulated data 
% to be used for modelling 

% ---------------

% remove trials of 50% prob (if there are any)!!
% ii50 = xstate==0.5; 
% xstate(ii50) = [];

% compute choice probability using softmax 
dv  = [0; val(1:end - 1)]-.5;
p   = 1./(1+exp(-dv));
% p(ii50) = [];

% define which are the stable and the volatile trials
trialstable     = tstable; 
trialvolatile   = tvolatile;

% need to work on this part as it seems that it is not working for our
% design
corr_action     = xstate>=.5;

choice          = p>=.5;
perf            = choice==corr_action;


mpvol = mean(perf(trialvolatile,:));
mpstab = mean(perf(trialstable,:));

mpvol = mean(mpvol);
mpstab = mean(mpstab);

dp = mpstab - mpvol;

xx = [mpstab mpvol];

end % end of function