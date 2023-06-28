function [xx, dp, trialstable, trialvolatile] = responseModel_v1(xstate, val,tvolatile, tstable)

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


% define which are the stable and the volatile trials
trialstable     = tstable; 
trialvolatile   = tvolatile;

% N = length(p);
% n = 10;
% change_points = find(diff(xstate)~=0);
% tvolatile = false(N,1);

% volatile trials are the trials with no probabilitiys witch in the preceding 10 trials
% for i=1:length(change_points)
%     tvolatile(change_points(i) + (1:n)) = 1;    
% end


% tstable = ~tvolatile;
% % tstable(1:n) = 0;


% need to work on this part as it seems that it is not working for our
% design
corr_action     = xstate>=.5;

choice          = p>=.5;
perf            = choice==corr_action;

mpvol           = mean(perf(trialvolatile,:));
mpstab          = mean(perf(trialstable,:));

mpvol           = mean(mpvol);
mpstab          = mean(mpstab);

dp              = mpstab - mpvol;
xx              = [mpstab mpvol];


end % end of function