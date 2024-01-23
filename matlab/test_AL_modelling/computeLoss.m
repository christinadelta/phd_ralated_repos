function [sum_simloss, sum_fitloss] = computeLoss(fitaction, simactions,simoutcome)

% this function computes loss for simulated and recovered loss

%% how many reps?

reps    = size(fitaction,2);
loss    = 0.5;
noloss  = 0;

%% extract actions and outcomes for each rep and compute loss

% loop over reps
for rep = 1:reps

    this_simaction  = simactions{1,rep};
    this_fitaction  = fitaction{1,rep};
    this_outcome    = simoutcome{1,rep};
    tmp_simloss     = 0; % init loss
    tmp_fitloss     = 0; % init loss

    % lets compute loss on a trial by trial basis 
    for i = 1:3 

        for j = 1:2

            tmp_simaction   = this_simaction{i,j};
            tmp_fitaction   = this_fitaction{i,j};
            tmp_outcome     = this_outcome{i,j};

            for k = 1:length(tmp_outcome)
                
                % first updated simulated loss
                if tmp_simaction(k) == tmp_outcome(k) % if good option is chosen
                    tmp_simloss(k) = noloss; % update loss
                else % if bad option is chosen
                    tmp_simloss(k) = loss; % update loss
                end

                % update recovered loss
                if tmp_fitaction(k) == tmp_outcome(k) % if good option is chosen
                    tmp_fitloss(k) = noloss; % update loss
                else % if bad option is chosen
                    tmp_fitloss(k) = loss; % update loss
                end  
            end % end of condition trials loop

            % store loss
            simloss{rep}{i,j} = tmp_simloss';
            fitloss{rep}{i,j} = tmp_fitloss';

        end % end of volatility loo
    end % end of stc loop
end % end of repetitions loop

%% sum loss for each repetition

for i = 1:reps 
    this_rep_simloss = simloss{1,i};
    this_rep_fitloss = fitloss{1,i};

    % re-arrange and sum simulated loss
    all_sim_loss(1,i) = sum(this_rep_simloss{1,1}); % small stc/stable
    all_sim_loss(2,i) = sum(this_rep_simloss{1,2}); % small stc/volatile
    all_sim_loss(3,i) = sum(this_rep_simloss{2,1}); % medium stc/stable
    all_sim_loss(4,i) = sum(this_rep_simloss{2,2}); % medium stc/volatile
    all_sim_loss(5,i) = sum(this_rep_simloss{3,1}); % large stc/stable
    all_sim_loss(6,i) = sum(this_rep_simloss{3,2}); % large stc/volatile

    % re-arrange and sum recovered loss
    all_fit_loss(1,i) = sum(this_rep_fitloss{1,1}); % small stc/stable
    all_fit_loss(2,i) = sum(this_rep_fitloss{1,2}); % small stc/volatile
    all_fit_loss(3,i) = sum(this_rep_fitloss{2,1}); % medium stc/stable
    all_fit_loss(4,i) = sum(this_rep_fitloss{2,2}); % medium stc/volatile
    all_fit_loss(5,i) = sum(this_rep_fitloss{3,1}); % large stc/stable
    all_fit_loss(6,i) = sum(this_rep_fitloss{3,2}); % large stc/volatile

    % sum loss over conditions
    sum_simloss(i) = sum(all_sim_loss(:,i));
    sum_fitloss(i) = sum(all_fit_loss(:,i));
end 


end % end of function