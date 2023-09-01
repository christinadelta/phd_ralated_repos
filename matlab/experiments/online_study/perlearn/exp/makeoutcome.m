function [tone_seq, outcome_seq] = makeoutcome(trlvec, prob)

% inputs: 
%           - trlvec = trials vector (e.g., 1:100)
%           - prob = high or low rewarded probability to be used to generate feedback sequence

% --------------------------------------

% create fixed feedback sequence 
cuetrls                         = length(trlvec)/2;
probtrls                        = round(cuetrls*prob);

% low tone outcome 
lowt_col(1:cuetrls,1)           = 1; % low tone
lowt_col(1:probtrls,2)          = 2; % house outcome
lowt_col(probtrls+1:cuetrls,2)  = 1; % face outcome

% high tone outcome 
hight_col(1:cuetrls,1)              = 2; % high tone
hight_col(1:probtrls,2)             = 1; % face outcome
hight_col(probtrls+1:cuetrls,2)     = 2; % house outcome

tone_seq                            = cat(1,lowt_col(:,1),hight_col(:,1));
outcome_seq                         = cat(1,lowt_col(:,2),hight_col(:,2));


end % end of function