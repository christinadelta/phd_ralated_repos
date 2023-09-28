function [tone, feed] = makeoutcome_v2(trlvec, prob, j)

% inputs: 
%           - trlvec = trials vector (e.g., 1:100)
%           - prob = high or low rewarded probability to be used to generate feedback sequence

% --------------------------------------

% create fixed cue/outcome sequence 
trialnum    = length(trlvec);
vertrials   = round(trialnum*prob);

% create cue/outcome list 
if j == 2                                   % if stc is medium the p(house|high tone)= high and p(face|low tone) = low
    tone(1:vertrials,1)             = 0;    % high tone - high prob
    tone(vertrials+1:trialnum,1)    = 1;    % low tone - low prob
    feed(1:vertrials,1)             = 1;    % house outcome - high probability 
    feed(vertrials+1:trialnum,1)    = 0;    % face outcome - low probability 

else                                        % if stc is small or large then, p(house|lowtone)= high and p(face|low tone) = low
    tone(1:vertrials,1)             = 1;    % low tone - high prob
    tone(vertrials+1:trialnum,1)    = 0;    % high tone - low prob
    feed(1:vertrials,1)             = 1;    % house outcome - high probability 
    feed(vertrials+1:trialnum,1)    = 0;    % face outcome - low probability 
end


end % end of function