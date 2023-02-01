function [performance] = heatmap_combparams(tmp_params,cond,probs,trials,outpath,task)

% the function is used to simulate, data run the RW model and average
% choice probabilities and performance across all combinations of alpha and
% beta parameter values

% --------------------------------------
% simulate data 
data        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output

% run model
modelout    = modelRW_v1(tmp_params, data, outpath);

% extract performance (% correct)
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
performance = nanmean(tmptrial);

% extract choice probabilities 
% choiceprobv = nanmean(modelout.allPs(:,1));
% choiceprobh = nanmean(modelout.allPs(:,2));


end % end of function