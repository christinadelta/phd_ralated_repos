function [trlbytrl_choices] = combparams(condition,probs,trials, outpath,task)
% function [trlbytrl_choices, correctchoicese, correctchoicesl] = combparams(condition,probs,trials, outpath,task)

% created 29/01/2023

% comments should be included here
% describe function


% ---------------------------------------------------
% use different parameter values and re-run the model 
alphas  = [0.20 0.5 0.75 1]; % use different alpha values and re-run the model 
betas   = [3 5 9 15];

for rep = 1:1000
    for alpha = 1:length(alphas)
    
        for beta = 1:length(betas)
    
            for cond = 1:condition
    
                temp_params = [alphas(alpha) betas(beta)]; 
                
                % simulate data 
                data        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
                
                % run model
                modelout    = modelRW_v1(temp_params, data, outpath);
                
                % extract feedback probabilities and actions to compute
                % correct choices 
                tmpchoice       = data.feedbackprob;
                tmpchoice(:,2)  = modelout.a';
                
                % extract actions/choices for the high probability
                if cond == 1 % if stable condition

                    if task == 2 % if stable with one switch
                        for i = 1:trials

                            if tmpchoice(i,1) == probs(1,1) && tmpchoice(i,2) == 1

                                tmptrial(i) = 1;
                            elseif tmpchoice(i,1) == probs(1,1) && tmpchoice(i,2) == 2
                                tmptrial(i) = 0;
                                
                            elseif tmpchoice(i,1) == probs(1,2) && tmpchoice(i,2) == 2

                                tmptrial(i) = 1;
                            else
                                tmptrial(i) = 0;
                            end
                        end
                    else

                        [~, imax] = max(probs(1,:));
                        tmptrial = modelout.a == imax;

                    end
                else % if condition = 2 (volatile condition

                    for i = 1:trials

                        if tmpchoice(i,1) == probs(2,1) && tmpchoice(i,2) == 1

                            tmptrial(i) = 1;
                        elseif tmpchoice(i,1) == probs(2,1) && tmpchoice(i,2) == 2
                            tmptrial(i) = 0;
                            
                        elseif tmpchoice(i,1) ~= probs(2,1) && tmpchoice(i,2) == 2

                            tmptrial(i) = 1;
                        else
                            tmptrial(i) = 0;
                        end
                    end
               
                end % end of if statement

                % store trial-by-trial vertical-gabor choices for ploting 
                trlbytrl_choices{cond}{alpha}{beta}(rep,:) = tmptrial;
                
                
                % store correct choices for the vertical gabor (high
                % prob gabor) for all combinations of alphas and betas
                % correctchoices{cond}(alpha,beta,rep) = nanmean(modelout.a == imax);
%                 correctchoicese{cond}(alpha,beta,rep) = nanmean(modelout.a(1:10) == imax); % only get the first 10 trials
%                 correctchoicesl{cond}(alpha,beta,rep) = nanmean(modelout.a(end-9:end) == imax); % only get the last 10 trials

            end % end of condition loop
    
        end % end of betas

    end % end of alphas loop
end % end of repetition loop


end % end of function
