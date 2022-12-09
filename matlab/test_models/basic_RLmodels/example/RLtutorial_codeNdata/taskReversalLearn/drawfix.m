% draw a fixation cross, consisting of 2 lines
%
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 02-08-2012
%==========================================================================

wdw = prep.draw.wdw;
wdh = prep.draw.wdh;
fx  = prep.draw.fx;
Screen('DrawLine',wd,prep.draw.grey,wdw/2, wdh/2-fx, wdw/2, wdh/2+fx,2);
Screen('DrawLine',wd,prep.draw.grey,wdw/2-fx, wdh/2, wdw/2+fx, wdh/2,2);
