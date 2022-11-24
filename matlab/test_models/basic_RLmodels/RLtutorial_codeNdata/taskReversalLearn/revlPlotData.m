function fh = revlPlotData(subjects,dataPath)
% FUNCTION revlPlotdata(subjects, dataPath)
% 
% Plot the raw behavioural data and smoothed response function for the the
% reversal learning task. 
% 
% subjects = subject ID vector, e.g. [1:2:10] 
% dataPath = full path to directoroy that contains the data
% 
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 11-08-2015
%==========================================================================


volatility  = {'high','low'};

for sID = subjects(:)'
    taskVersion = 2-mod(sID,2); % determine which version of the task they did
    
    % define paths and filenames
    subTag      = sprintf('tutorialRevLearn_%s_s%03.0f',volatility{taskVersion},sID);
    dataFile = fullfile(dataPath,sprintf('%s_data.mat', subTag));
    load(dataFile) % get data
    
    % recode choices so they're 0/1 where 1 is the initially correct
    % stimulus
    choice = 2-data.choice(:)'; 
    smoothingkernel = 6; % average across 7 trials
    
    % plot everything
    fh = figure; box off; hold on;
    set(fh,'position', [10 60 700 400],'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
    title(sprintf('subject %02.0f, %s volatility',sID,volatility{taskVersion}));
    plot(data.prep.feedbackprob,'k-','linewidth',2);
    plot(choice,'k*','markersize',10);
    plot(mySmooth(choice,smoothingkernel,[],'backward'),'r:','linewidth',2)
    ylabel('p(reward|blue)');
    ylim([-0.1 1.1]);
    xlabel('trial');
    legend({'p(reward|blue)','choice','smoothed choice'},'location','northeastoutside')
    legend boxoff
end

