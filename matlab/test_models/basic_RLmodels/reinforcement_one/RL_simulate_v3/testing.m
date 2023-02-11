if condition == 1 && task == 1
    tmp_feedbackv       = cat(1,ones(runtrials * prob,1), zeros(runtrials * (1-prob),1));
    shuffled            = randperm(length(tmp_feedbackv));% shuffle the feedback 
    feedbackv           = tmp_feedbackv(shuffled);
    
    clear shuffled tmp_feedbackv
    % let's create an independent column for the second stimulus (horizontal
    % gabor) 
    tmp_feedbackh       = cat(1,zeros(runtrials * prob,1), ones(runtrials * (1-prob),1));
    shuffled            = randperm(length(tmp_feedbackh));% shuffle the feedback 
    feedbackh           = tmp_feedbackh(shuffled);
    
    % add columns two one matrix
    % feedback             = cat(2, feedbackv, feedbackh);

elseif condition == 1 && task == 2

    % generate column 1 (vertical gabor reward feedback)
    % first create the 75% trials (for vertical gabor)
    tmp_feedbackv1       = cat(1,ones(ceil(runtrials * prob),1), zeros(floor(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackv1));% shuffle the feedback 
    feedbackv1           = tmp_feedbackv1(shuffled1);

    % now create the 25% trials (for vertical gabor)
    tmp_feedbackv2       = cat(1,ones(ceil(runtrials * (1-prob)),1), zeros(floor(runtrials * prob),1));
    shuffled2            = randperm(length(tmp_feedbackv2));% shuffle the feedback 
    feedbackv2           = tmp_feedbackv2(shuffled2);

    % concatinate the two arrays
    feedbackv            = cat(1, feedbackv1, feedbackv2);

    clear shuffled1 shuffled2

    % generate column 2 (horizontal gabor reward feedback)
    % first create the 75% trials (for horizontal gabor)
    tmp_feedbackh1       = cat(1,zeros(ceil(runtrials * prob),1), ones(floor(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackh1));% shuffle the feedback 
    feedbackh1           = tmp_feedbackh1(shuffled1);

    % now create the 25% trials (for horizontal gabor)
    tmp_feedbackh2       = cat(1,ones(ceil(runtrials * prob),1), zeros(floor(runtrials * (1-prob)),1));
    shuffled2            = randperm(length(tmp_feedbackh2));% shuffle the feedback 
    feedbackh2           = tmp_feedbackh2(shuffled2);

    % concatinate the two arrays
    feedbackh            = cat(1, feedbackh1, feedbackh2);

    clear shuffled1 shuffled2

else % if this is volatile condition

    % first generate feedback for column 1 (vertical gabor)
    % 80% feedback probability
    tmp_feedbackv1      = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled1           = randperm(length(tmp_feedbackv1)); % shuffle the feedback
    feedbackv1          = tmp_feedbackv1(shuffled1);

    % 20% feedback probability
    tmp_feedbackv2      = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled2           = randperm(length(tmp_feedbackv2)); % shuffle the feedback
    feedbackv2          = tmp_feedbackv2(shuffled2);

    % 80% feedback probability
    tmp_feedbackv3      = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled3           = randperm(length(tmp_feedbackv3)); % shuffle the feedback
    feedbackv3          = tmp_feedbackv3(shuffled3);

    % 20% feedback probability
    tmp_feedbackv4      = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled4           = randperm(length(tmp_feedbackv4)); % shuffle the feedback
    feedbackv4          = tmp_feedbackv4(shuffled4);

    % concatinate the two arrays
    feedbackv            = cat(1, feedbackv1, feedbackv2, feedbackv3, feedbackv4);

    clear shuffled1 shuffled2 shuffled3 shuffled4

    % now generate feedback for column 2 (horizontal gabor)
    % 80% feedback probability
    tmp_feedbackh1       = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackh1));% shuffle the feedback 
    feedbackh1           = tmp_feedbackh1(shuffled1);

    % 20% feedback probability
    tmp_feedbackh2       = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled2            = randperm(length(tmp_feedbackh2));% shuffle the feedback 
    feedbackh2           = tmp_feedbackh2(shuffled2);
    
    % 80% feedback probability
    tmp_feedbackh3       = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled3            = randperm(length(tmp_feedbackh3));% shuffle the feedback 
    feedbackh3           = tmp_feedbackh3(shuffled3);

    % 20% feedback probability
    tmp_feedbackh4       = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled4            = randperm(length(tmp_feedbackh4));% shuffle the feedback 
    feedbackh4           = tmp_feedbackh4(shuffled4);

    % concatinate the two arrays
    feedbackh            = cat(1, feedbackh1, feedbackh2, feedbackh3, feedbackh4);

    clear shuffled1 shuffled2 shuffled3 shuffled4

end

% add columns to one matrix
feedback             = cat(2, feedbackv, feedbackh);





%% work on Nura's code:

ind = 1:100;
prob = 0.75;
chance = 1;
rewVals = [1 -1];
