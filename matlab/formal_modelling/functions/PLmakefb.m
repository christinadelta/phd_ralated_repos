function [tone_seq, outcome_seq] = PLmakefb(trlvec, prob)

% inputs: 
%           - trlvec = trials vector (e.g., 1:100)
%           - prob = high or low rewarded probability to be used to generate feedback sequence

% --------------------------------------

% create fixed feedback sequence 
cuetrls                         = length(trlvec)/2;
probtrls                        = round(cuetrls*prob); % just to be safe

% low tone outcome 
lowt_col(1:cuetrls,1)           = 1; % low tone
lowt_col(1:probtrls,2)          = 0; % house outcome
lowt_col(probtrls+1:cuetrls,2)  = 1; % face outcome

% high tone outcome 
hight_col(1:cuetrls,1)              = 0; % high tone
hight_col(1:probtrls,2)             = 1; % face outcome
hight_col(probtrls+1:cuetrls,2)     = 0; % house outcome

tone_tmp                            = cat(1,lowt_col(:,1),hight_col(:,1));
outcome_tmp                         = cat(1,lowt_col(:,2),hight_col(:,2));
seq(:,1) = tone_tmp;
seq(:,2) = outcome_tmp;

% shuffle 
randomiser  = randperm(length(tone_tmp));
seq         = seq((randomiser),:);

tone_seq = seq(:,1);
outcome_seq = seq(:,2);



end % end of function