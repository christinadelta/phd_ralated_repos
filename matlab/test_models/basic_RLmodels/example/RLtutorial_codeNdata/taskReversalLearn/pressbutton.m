%==========================================================================
% PRESSBUTTON
% script to display and wait for a keypress to continue
%
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 11-08-2015
%==========================================================================

tekst = 'Press any key to continue';
[wt]=Screen(wd,'TextBounds',tekst);
Screen('Drawtext',wd,tekst,wdw*0.38,wdh-(prep.top+2*prep.draw.txtsize),white);
Screen('Flip', wd);

WaitSecs(0.2)
while true % using keyboard.
    [~,keyCode,~]= KbWait;
    if (find(keyCode) == prep.abort) % abort
        Screen('CloseAll');
        break;
    else
        break
    end
end
