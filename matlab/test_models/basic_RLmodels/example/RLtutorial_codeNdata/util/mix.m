%==========================================================================
% Returns a randomised version of the input vector. 
%
% Hanneke den Ouden,
% start:        16-08-2009
% last changes  02-08-2012
%==========================================================================

function mix = mix(A)
mix = A(randperm(length(A)));
return