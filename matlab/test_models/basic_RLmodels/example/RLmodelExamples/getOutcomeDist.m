function [out] = getOutcomeDist(ind, prob, chance, rewVals)
% Create the distribution of positive & negative outcomes
% input--   ind: indexes for relevant condition
%           prob: reward probability to apply
%           chance: if 1 (default), create reward distribution stochastically (given prob & rng)
%                   if 0, distribution is fixed based on n
%           rewVals: values for rewarding and non-rewarding outcomes, default=[1, -1] 
%%% NS Apr 2020
    
if nargin < 4 || isempty(rewVals)
    rewVals = [1, -1];
end
if nargin < 3 || isempty(chance)
    chance = 1;
end

if numel(ind) == 1 % input required number of trials
    n = ind;
else               % input is vector of trial indices 
    n = length(find(ind));
end


shufflef = @(v)v(randperm(numel(v)));

out = NaN(n, 1);

if ~chance    
    nGood = prob*n;
    if ~rem(nGood,1) % it's a whole number
        out(1:nGood) = rewVals(1); % Good reward
        out(nGood+1:n) = rewVals(2); % Bad reward       
        out = shufflef(out); % randomise order                        
        
    else    %%% If prob*n is not a whole number, probabilities would be
            %%% systematically biased by the rounding, hence,
            %%% must set to stochastic draw or change n
        error(' getOutcomeDist: inconsistent input -- prob=%4.4f * nTrialCd=%d does not yield integer, but chance=0 - i.e. calls for non-stochastic assignment.', prob, n)
    end
    
else % stochastic draw, so order is already "randomised"
    for i = 1:n
        if rand <= prob
            out(i) = rewVals(1); % Good reward
        else
            out(i) = rewVals(2); % Bad reward
        end
    end    
end