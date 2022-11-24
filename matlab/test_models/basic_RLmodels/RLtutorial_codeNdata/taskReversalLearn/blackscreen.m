% blackscreen
% Bit of script to show a black screen in between two images
% 
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 02-08-2012
%==========================================================================

Screen('FillRect',wd,black);
Screen('Flip', wd);
WaitSecs(0.05);
