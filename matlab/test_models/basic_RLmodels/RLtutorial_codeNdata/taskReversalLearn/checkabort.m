% checkabort
% check whether abort key was pressed
% 
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 02-08-2012
%==========================================================================

[keyIsDown,~, keyCode]= KbCheck;
if keyIsDown && any(find(keyCode) == prep.abort)
    aborted = true;
    break
end
clear keyIsDown secs keyCode

