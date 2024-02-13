function [meanNLL,allrepsNLL] = getMeanNLL(allNLL)

%% HOW MANY REPETITIONS?
nreps = size(allNLL,2);

for rep = 1:nreps

    tempNLL = allNLL{1,rep};

    allrepsNLL(rep,1) = tempNLL(1,1); % NLL small/stable
    allrepsNLL(rep,2) = tempNLL(1,2); % NLL small/volatile
    allrepsNLL(rep,3) = tempNLL(2,1); % NLL medium/stable
    allrepsNLL(rep,4) = tempNLL(2,2); % NLL medium/volatile
    allrepsNLL(rep,5) = tempNLL(3,1); % NLL large/stable
    allrepsNLL(rep,6) = tempNLL(3,2); % NLL large/volatile

    clear tempNLL
end % end of reps loop

%% average nll

meanNLL = mean(allrepsNLL,1);


end % end of function