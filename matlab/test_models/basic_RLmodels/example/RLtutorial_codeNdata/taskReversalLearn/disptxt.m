function time = disptxt(txt,wd,wdw,wdh,align, top, col,wait,dontblank)
% FUNCTION time = disptxt(txt,wd,wdw,wdh,align, top, col,wait,dontblank)
% Display text, with various alignments determined by the input parameters
% 
% txt = text to be displayed
% wd = window pointer
% wdw = window width
% wdh = window height
% align = left-align text if true, centre otherwise
% top = top-align text if true, centre otherwise
% col = text colour
% wait = wait for .2 seconds before subjects can respond if true
% dontblank = if true, blank before flipping
% 
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 11-08-2015
%==========================================================================

nrow=length(txt);
for k=1:nrow
    [wt] = Screen('TextBounds',wd,txt{k});
    if align
        xpos = wdw*0.2;
    else %centre
        xpos = round(wdw/2-wt(3)/2);
    end
    if top
        ypos = round(wdh/12+2*(k-1)*wt(4));
    else %centre
        ypos = round(wdh/2+2*(k-1-nrow/2)*wt(4));
    end
    Screen('Drawtext',wd,txt{k},xpos,ypos,col);
end

time = GetSecs;
if dontblank; Screen('flip', wd,[],1);
else         Screen('flip', wd);
end

if wait;  WaitSecs(0.2); KbWait;end
