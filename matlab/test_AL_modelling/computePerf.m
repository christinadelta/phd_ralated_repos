function modelPerf = computePerf(temp_outcome, temp_action)

% compute performance of fitted choices using participant outcomes 



% loop over stochasticity and volatility levels
for i = 1:3
    for j = 1:2


        % extract outcomes for good options and actions
        this_outcome                        = temp_outcome{i,j}(:,1);
        this_outcome(find(this_outcome==0)) = 2; % convert 0s to 2s 
        this_choice                         = temp_action{i,j};

        % how many reps?
        reps = size(this_choice,2);

        for k = 1:reps

            action  = this_choice(:,k);
            correct(k) = mean(this_outcome == action);

        end % end of reps loop

        all_correct = mean(correct);
        modelPerf(i,j)   = all_correct;

    end % end of volatility loop
end % end of stc loop


%%

end % end of function