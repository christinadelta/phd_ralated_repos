function [trlbytrl_choices] = combparams(condition,probs,trials, outpath,task)
% function [trlbytrl_choices, correctchoicese, correctchoicesl] = combparams(condition,probs,trials, outpath,task)

% created 29/01/2023

% comments should be included here
% describe function


% ---------------------------------------------------
% use different parameter values and re-run the model 
alphas  = [0.25 0.5 0.75 1];
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
                
                if cond == 1
                    if task == 1
                        [~,imax]    = max(probs(1,:));
                        tmptrial    = modelout.a == imax;
                    end
                else % if condition == 2
                    runs            = 4; 
                    trlruns         = trials/runs;
                    counter         = 1;

                    for run = 1:runs
                        
                        % work with p(correct) for the volatile condition
                        tmp                     = modelout.a(counter:trlruns*run);

                        if mod(run,2)
                            [~,imax]            = max(probs(2,:));
                            trls{run}           = tmp == imax;
                        else
                            [~,imin]            = min(probs(2,:));
                            trls{run}           = tmp == imin;
                        end

                        counter                 = counter + trlruns;
                    end % end of for loop
                    
                    tmptrial                    = [trls{1} trls{2} trls{3} trls{4}]; % concatenate

                end

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
