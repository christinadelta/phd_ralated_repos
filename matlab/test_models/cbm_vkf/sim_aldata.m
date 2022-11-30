
% create vectors with data needed for now..
% will keep building on this 

% data needed to try simplied versions of models:
% actions 
% outcome (correct vs incorrect)
% create a vector with actions (for now consider only face and house without tone association) just for the start 
stimcat = repmat([1 2], 23, 2); % create a list with ones=face and twos=house
stimcat = stimcat(:); % 
stimcat = stimcat(randperm(length(stimcat))); % randomise the sequence of stimuli
actions = stimcat;


% create vector with outcomes 
% ones =  



