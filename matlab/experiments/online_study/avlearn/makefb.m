function feedback = makefb(trlvec, prob, rdm)

% inputs: 
%           - trlvec = trials vector (e.g., 1:100)
%           - prob = high or low rewarded probability to be used to generate feedback sequence
%           - rdm = [if 1 = feedback sequence generation is stochastic (like all tutorials), if 0 = sequence is fixed based on the number of trials for high and low rewarded trials]      
%           - volatility [1x1 double, needed for the generation of the feedback sequence
%           - task (if 1 = in stable condition there is a switch of reward probabilities and the stimuli after half trials) 

% --------------------------------------
shuffleseq  = @(v)v(randperm(numel(v))); % for shuffling the sequence
trialnum    = length(trlvec);

% if rdm = 0, create fixed feedback distribution 
if rdm == 0 % if generate fixed sequence
    vertrials                           = round(trialnum*prob);
    feedback(1:vertrials,1)             = 1; % high probaility option
    feedback(vertrials+1:trialnum,1)    = 2; % low probability option
    feedback                            = shuffleseq(feedback);

else % if generate random sequence
    for trl = 1:trialnum
        if double(rand(1) <= prob)
            feedback(trl,1)             = 1;
        else
            feedback(trl,1)             = 2;
        end
    end % end of for loop
end

end % end of function