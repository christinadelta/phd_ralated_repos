function yy = mySmooth(varargin)
% FUNCTION yy = mySmooth(varargin)
% 
% Smooth response data across the second dimension
% Syntax
% yy = smooth(y)
% yy = smooth(y,span)
% yy = smooth(y,span,dim)
% 
% Description
% 
% yy = smooth(y) smooths the data in the column vector y using a moving
% average filter. Results are returned in the column vector yy. 
% varargin{2} = span for moving average. The default is 5.  
% varargin{3} = dimension across which to smooth
% varargin{4} = how to smooth. Options are: 'backward' or 'centre':
% smoothing kernel centred on the current timepoint, or only smoothing back
% in time. Default = 'centre'; 
% 
% The first few elements of yy are given by
% 
% yy(1) = y(1)
% yy(2) = (y(1) + y(2) + y(3))/3
% yy(3) = (y(1) + y(2) + y(3) + y(4) + y(5))/5
% yy(4) = (y(2) + y(3) + y(4) + y(5) + y(6))/5
% 
% yy = smooth(y,span) sets the span of the moving average to span. span
% must be odd. 
% 
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2012-2015 <h.denouden@donders.ru.nl>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------

y = varargin{1};
if numel(varargin)<2
    span = 5;
elseif isempty(varargin{2})
    span = 5;
else 
    span = varargin{2};
end

if numel(varargin)<3
    dim = 2;
elseif isempty(varargin{3})
    dim = 2;
else
    dim = varargin{3};
end

if numel(varargin)<4
    direction = 'center';    
else
    direction = varargin{4};    
end

if dim==1
    y = y';
end

% check if span is odd
if ~mod(span,2) && strcmp(direction,'center')
    error('if using centred smoothing, span needs to be an odd number!')    
end

% check whether the 
nt = length(y);
if ~any(size(y)==1)
[ns, nt] = size(y);
end

yy = nan(size(y));

if strcmp(direction,'center')
    k = floor(span/2);
    for t = 1:nt
        if t<=k
            c = floor(t/2);
        elseif t>k && t<nt-k
            c = k;
        else
            c = nt-t;
        end
        yy(:,t) = mean(y(:,t-c:t+c),2);
    end
elseif strcmp(direction,'backward')
    k = span;
    for t = 1:nt
        if t<k
            c = t-1;
        else
            c = k-1;
        end
        yy(:,t) = mean(y(:,t-c:t),2);
    end
end


if dim ==1
    yy = yy';
end

