function [stimuli_left, stimuli_right] = writestimuli(blockTrials)

% create cell of stimuli for stable and volatile blocks

% stable stimuli
stableStim = {'blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png'...
    'blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png'...
    'red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png'...
    'red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png'}';

% volatile stimuli
volStim = {'blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png'...
    'blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png'...
    'blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png','blue.png'...
    'red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png'...
    'red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png'...
    'red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png','red.png'}';


% now for each block copy the celss and randosize
%v = [v1,v2];
shufflestim  = @(v)v(randperm(numel(v)));
one = shufflestim(stableStim);
two = shufflestim(volStim);
three = shufflestim(stableStim);
four = shufflestim(stableStim);
five = shufflestim(volStim);
six = shufflestim(stableStim);
seven = shufflestim(stableStim);
eight = shufflestim(volStim);
nine = shufflestim(stableStim);


%% concatinate them all in one array of stimuli

stimuli_left = [one;two;three;four;five;six;seven;eight;nine];

%% make the second array with the opposite stimuli

% if stim left is 'blue' stim right should be 'red'
l = length(stimuli_left);

for i = 1:l

    if length(stimuli_left{i}) == 8 % if blue 

        stimuli_right{i} = 'red.png';
    elseif length(stimuli_left{i}) == 7 % if red
        stimuli_right{i} = 'blue.png';
    end
end

stimuli_right = stimuli_right';

%% concatinate the left and right arrays
% allstimuli = [stimuli_left, stimuli_right]

end

