function [stable_corr, vol_corr] = getCorrect(o,actions)

% get performance for model simulated responses 

% loop over stc levels
for i = 1:3
    for j = 1:2

        action_tmp  = actions{i,j};
        o_tmp       = o{i,j}(:,1,:);

        for sim = 1:size(action_tmp,2)

            a_sim = action_tmp(:,sim);
            o_sim = o_tmp(:,1,sim);

            % loop over trials 
            for t = 1:length(o_sim)

                if a_sim(t,1) == o_sim(t,1)
                    correct(t,sim) = 1;
                else
                    correct(t,sim) = 0;
                end
            end % end of trials loop
        end % end of simulations loop

        % average performance across simulations 
        temp_correct = mean(correct,1)';

        % add performance to the correct array
        if j == 1
            stable_corr(:,i) = temp_correct;
        elseif j == 2
            vol_corr(:,i) = temp_correct;
        end

    end % end of vol loop
end % end of stc loop

end % end of function 