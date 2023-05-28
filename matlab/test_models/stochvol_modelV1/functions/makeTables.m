function [lowStochTable highStochTable] = makeTables(x,y)

% make tables to visualise the simulated data

%% low stochasticity table

% loop over volatility condition
for i = 1:2

    StochProbRel{1,i} = y{1,i};
    outcome{1,i}(:,1) = x{1,1}{i,1};
    outcome{1,i}(:,2) = x{1,1}{i,2};

end % end of volatility loop

% concatenate the arrays
lowProbRel = [StochProbRel{1,1}; StochProbRel{1,2}];
lowOutcome = [outcome{1,1}; outcome{1,2}];

% make table
lowStochTable = table(lowProbRel, lowOutcome);

clear outcome StochProbRel

%% high stochasticity table

% loop over probability condition
for i = 1:2

    StochProbRel{1,i}   = y{2,i};
    outcome{1,i}(:,1)   = x{1,2}{i,1};
    outcome{1,i}(:,2)   = x{1,2}{i,2};

end % end of volatility loop

highProbRel             = [StochProbRel{1,1}; StochProbRel{1,2}];
highOutcome             = [outcome{1,1}; outcome{1,2}];
highStochTable          = table(highProbRel, highOutcome);


end % end of function