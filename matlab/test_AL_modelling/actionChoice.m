function a = actionChoice(p)

% to be used used in the response model to chose an action based on the
% probabilities computed with the softmax function
% ---------

a = max(find([-eps cumsum(p)] < rand));

end  % end of function

