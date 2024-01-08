function [feedback, cue] = makefb_v2(trlvec, prob, rdm)

% inputs: 
%           - trlvec = trials vector (e.g., 1:100)
%           - prob = high or low rewarded probability to be used to generate feedback sequence
%           - rdm = [if 1 = feedback sequence generation is stochastic (like all tutorials), if 0 = sequence is fixed based on the number of trials for high and low rewarded trials]      

% --------------------------------------
trialnum    = length(trlvec);
shuffleseq  = @(v)v(randperm(numel(v))); % for shuffling the sequence

% if rdm = 0, create fixed feedback distribution 
if rdm == 0 % if generate fixed sequence
    vertrials                           = round(trialnum*prob);
    feedback(1:vertrials,1)             = 1; % high probaility option
    feedback(vertrials+1:trialnum,1)    = 0; % low probability option
    cue(1:vertrials,1)                  = 1; % high probaility option
    cue(vertrials+1:trialnum,1)         = 0; % low probability option
    feedback                            = shuffleseq(feedback);


else % if rdm == 1, generate random sequence
    for trl = 1:trialnum
        if double(rand(1) <= prob)
            feedback(trl,1)             = 1;
            cue(trl,1)                  = 1;
        else
            feedback(trl,1)             = 0;
            cue(trl,1)                  = 0;
        end
    end % end of for loop
end

end % end of function