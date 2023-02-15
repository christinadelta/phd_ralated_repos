function fitted = RLtutorial_main
% ------------------------------------------------------------------------
% FUNCTION fitted = RLtutorial_main
%
% Tutorial to do simple modelling, associated with the online tutorial and
% dataset for a simple reversal learning tasks with 2 levels of volatility.
%
% Using this code, you can factorially vary whether to
% - simulate data or use real data
% - just visualise the (simulated or real) data or estimate the model
% - plot individual data or data aggregated over subjects
%
% The learning model used is a simple Rescorla-Wagner (Rescorla & Wagner
% 1972) function, using a softmax choice/link/observation function:
%   learning model (RW):     vA <- vA + alpha*(r-vA)
%   choice function (softmax):    p(A) = exp(beta*vA)/(exp(beta*vA)+exp(beta*vB))
%
% The parameter estimation is performed using a grid search. The function
% returns both expected value and maximum likelihood paramer estimates.
%
% Any parameters in the 'modify' section can be altered.
%
%
% Code by Hanneke den Ouden <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
%
% version: 11 August 2015
% ------------------------------------------------------------------------


%==========================================================================
%% Section 1: Preparation
%==========================================================================
% close all;
clear all
clc

%%%%%%%%%%%%    MODIFY      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjects        = 'all'; % use vector of numbers if simulating or using a subset of the data; use string 'all' if using all available real datasets
simulate        = false;    % true if simulate, false if using real data
fitData         = true; 	% if false, we only simulate / plot data
plotIndividual  = false;  % only set to to true if running a small number of individuals!

% If simulating data, specify the parameters you want to use for fitting
if simulate
    simpars.alpha = .25;
    simpars.beta  = 4;
end

% If fitting data, specify the bounds of the grid, and the number of bins
% that we want to compute.
if fitData
    bounds      = [0 1; % alpha
        0 15]; % beta
    nBin        = [20 25] ;
end

%%%%%%%%%%%%    MODIFY      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set paths for code we need
tmp = fileparts(which('RLtutorial_main'));
rootdir = tmp(1:end-11); % directory in which you put RL_tutorial_main
addpath(fullfile(rootdir,'ReinfLearn')); % modelling code
addpath(fullfile(rootdir,'util')); % some helper functions
addpath(fullfile(rootdir,'taskReversalLearn')); % to get the feedback code
dataPath = fullfile(rootdir,'taskReversalLearn','data'); % where our data comes from

% set labels that we'll need later on.
volatility  = {'high','low'};
paramLabel  = {'alpha','beta'};
nParam      = length(paramLabel);

% if using all subjects
if strcmp(subjects,'all')
    subFind = dir(fullfile(dataPath,'*_data*')); % select proper files
    subFind = {subFind.name};
    for iSub = 1:numel(subFind)
        idx = regexp(subFind{iSub},'\d'); % recognise numbers.
        if idx>0
            subNr = subFind{iSub}(idx); % select the number
            subNr = str2double(subNr);
            subList(iSub) = subNr;
        end
    end
    subjects = sort(subList); % create a final subjectlist
end

% sort subjects by task condition (high or low volatility);
highVol = mod(subjects,2)==1;
sortedsubs{1} = subjects(highVol); % odd
sortedsubs{2} = subjects(~highVol); % even
nsubTask = max(sum(highVol),length(subjects)-sum(highVol));
nTrial = 135;

%==================================================================
%% Section 2: Simulate data or plot real data
%==================================================================

% initialise variables
score           = nan(2,nsubTask);
if ~plotIndividual; fhData = figure('Name','Data');
    set(fhData,'position',[10 60 900 650],'paperunits','centimeters','Color','w');
else
    trialfh = nan(2,nsubTask); % figure handle for individual trialwise plots
end

for task = 1:2 % we'll do the subjects separately by task condition
    ct      = 0; % reset subject counter
    choice  = nan(nsubTask,nTrial);
    VV      = nan(nsubTask,nTrial,2);
    PP      = nan(nsubTask,nTrial);
    
    for sID = sortedsubs{task}
        ct = ct+1; % ubdate for each subject in task condition
        
        % Simulate data
        %==================================================================
        if simulate
            % Simulate & store data
            [data pout{task,ct}] = RLtutorial_simulate([simpars.alpha simpars.beta],sID);
            alldata{task,ct} = data;
            VV(ct,:,:) = pout{task,ct}.VV;
            PP(ct,:) = pout{task,ct}.PP(:,1);
            
            % Everything below if for plotting only
            choice(ct,:) = 2-data.choice; % convert to 0/1 for plotting
            if plotIndividual
                trialfh(task,ct) = figure; box off; hold on;
                set(trialfh(task,ct),'position', [10 60 700 400],'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
                plot(data.prep.feedbackprob,'k','linewidth',2)
                plot(pout{task,ct}.PP(:,1),'b-','linewidth',2)
                plot(pout{task,ct}.VV(:,1),':','color',[.4 0.4 1],'linewidth',2)
                plot(pout{task,ct}.VV(:,2),':','color',[1 .35 0.1],'linewidth',2)
                plot(2-pout{task,ct}.data.choice,'k*')
                legend({'p(reward|blue)','p(choose blue)','value(blue)','value(orange)',...
                    'choice'},'location','northeastoutside');
                legend boxoff
                ylabel('probability');
                ylim([-0.1 1.1]);
                xlabel('trial');
                title(sprintf('simulated data, alpha = %0.02f, beta = %02.01f',...
                    simpars.alpha,simpars.beta));
            end
            
            % ... or load and plot real data
            %==============================================================
        else
            subTag      = sprintf('tutorialRevLearn_%s_s%03.0f',volatility{task},sID);
            dataFile    = fullfile(dataPath,sprintf('%s_data.mat', subTag));
            load(dataFile,'data');
            alldata{task,ct} = data;
            if plotIndividual
                trialfh(task,ct)= revlPlotData(sID,dataPath);
            end
            choice(ct,1:data.prep.nt) = 2-data.choice; % convert to 0/1 for plotting
        end
        score(task,ct) = (sum(choice(ct,data.prep.feedbackprob==data.prep.prob)==1)...
            +sum(choice(ct,data.prep.feedbackprob==(1-data.prep.prob))==0))/data.prep.nt;
    end % individual subject
    
    % plot the data averaged across all subjects, for each task separately if not plotting
    % individuals
    %==============================================================
    if ~plotIndividual
        figure(fhData)
        subplot(2,1,task);hold on
        plot(data.prep.feedbackprob,'k:','linewidth',2)
        plot(nanmean(choice,1),'k-','linewidth',1)
        if simulate
            plot(mean(PP,1),'b-','linewidth',2)
            plot(mean(squeeze(VV(:,:,1)),1),':','color',[.4 0.4 1],'linewidth',2)
            plot(mean(squeeze(VV(:,:,2)),1),':','color',[1 .35 0.1],'linewidth',2)
            legend({'p(reward|blue)','mean choice (blue)','p(choose blue)','value(blue)','value(orange)',...
                },'location','northeastoutside');
            title(sprintf('MEAN data, %s volatility, alpha = %0.02f, beta = %02.01f',...
                volatility{task},simpars.alpha,simpars.beta));
        else
            legend({'p(reward|blue)','mean choice (blue)'},'location','northeastoutside');
            title(sprintf('MEAN data, %s volatility',volatility{task}));
        end
        legend boxoff
        ylabel('probability');
        ylim([0 1]);
        xlabel('trial');
        ylabel('probability');
    end % if plotting group
end % task

% plot the scores
if ~plotIndividual
    fhScore = figure('Name','score');set(fhScore,'position',[10 60 250 400],'paperunits','centimeters',...
        'paperposition',[0 0 6 6],'Color','w');
    barScatter(score,[],volatility,true);
    set(gca,'xtick',1:2,'xticklabel',volatility);
    xlabel('volatility');
    ylabel('p(correct)');
    title('correct choice')
    hline(.5,'w:');
    ylim([0 1])
    box off
end


%==========================================================================
%% Section 3.
% Do a grid search to compute the likelihood functions of the
% parameters, within prespecified bounds
%==========================================================================
if fitData
    if size(bounds,2)~=nParam
        error('number of bounds and parameters don''t match')
    end
    
    % Initialise the output matrices for the fitted parameters.
    fitted.alpha.pdf    = nan(2,nsubTask,nBin(1));
    fitted.beta.pdf     = nan(2,nsubTask,nBin(2));
    fitted.alpha.ml     = nan(2,nsubTask);
    fitted.beta.ml      = nan(2,nsubTask);
    fitted.alpha.ev     = nan(2,nsubTask);
    fitted.beta.ev      = nan(2,nsubTask);
    
    % open some figures
    fhParamPDF = figure('Name','parameter PDF');
    set(fhParamPDF,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    if plotIndividual
        fhGrid{1} = figure('Name','Likelihoods high vol');
        set(fhGrid{1},'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        fhGrid{2} = figure('Name','Likelihoods low vol');
        set(fhGrid{2},'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        
    else
        fhParBar = figure('Name','Parameter estimates');
        set(fhParBar,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    end
    
    % Do the actual grid search
    % =========================================================================
    if simulate % compute where the stimulated alpha and beta are in binspace
        simAlpha_bin = simpars.alpha/bounds(1,2)*nBin(1);
        simBeta_bin = simpars.beta/bounds(2,2)*nBin(2);
    end
    
    ctall = 0; % subject counter
    for task = 1:2
        ct = 0; % subject x task counter
        for sID = sortedsubs{task}
            ct = ct+1;
            ctall = ctall+1;
            sID_string{task}{ct} = sprintf('sub %03.0f',sID);
            data = alldata{task,ct};
            
            for iParam = 1:nParam
                range = linspace(bounds(iParam,1),bounds(iParam,2),nBin(iParam)+1);
                p{iParam} = range(2:end); % stay just off the zero bounds
            end
            
            params = nan(1,2);
            for t = 1:nBin(1)
                params(1) = p{1}(t);
                for tt = 1:nBin(2)
                    params(2) = p{2}(tt);
                    [loglik(t,tt), foo]= RLtutorial_fitmodel(params,data);
                end
            end
            
            loglik = loglik-min(loglik(:)); % remove the minimum;
            lik = exp(loglik); % compute the likelihood, (rather than the log)
            
            % compute the marginal likelihoods for each parameter: for
            % alpha, sum the likelihoods across all bins of beta, and  for
            % beta, sum the likelihoods across all bins of alpha. Then to
            % compute the probability density, normalise by dividing by the sum
            % of the likelihoods for each parameter.
            for x = 1:length(params)
                tmp = sum(lik,3-x);
                marglik{x} = tmp/sum(tmp);
            end
            
            % plot the likelihood landschape with the maximum and true
            % value (latter for simulations only), for each individual
            % =============================================================
            
            ML(1) = myvect(p{1}(max(marglik{1})==marglik{1}));
            ML(2) = myvect(p{2}(max(marglik{2})==marglik{2}));
            EV(1) = sum(p{1}(:).*marglik{1}(:));
            EV(2) = sum(p{2}(:).*marglik{2}(:));
            [foo, poutEst{task,ct,1}]= RLtutorial_fitmodel(ML,data);
            [foo, poutEst{task,ct,2}]= RLtutorial_fitmodel(EV,data);
            
            if plotIndividual
                % plot the grid
                figure(fhGrid{task});
                dims = ceil(sqrt(length(sortedsubs{task})));
                subplot(dims, dims,ct);
                imagesc(lik); hold on;
                
                [maxalpha maxbeta] = find(lik==max(lik(:)));
                hline(maxalpha, {'w-', 'LineWidth', 2});
                vline(maxbeta, {'w-', 'LineWidth', 2});
                
                if simulate
                    plot(simBeta_bin*ones(1,2),[simAlpha_bin-.75 simAlpha_bin+.75],'-','color',[0.4 0.4 0.4],'linewidth',3);
                    plot([simBeta_bin-.75 simBeta_bin+.75],simAlpha_bin*ones(1,2),'-','color',[0.4 0.4 0.4],'linewidth',3);
                end
                hold off
                
                title(sprintf('p(data|model), sub %03.0f, %s vol.',sID,volatility{task}));
                ytick = 0:nBin(1)/5:nBin(1);
                xtick = 0:nBin(2)/5:nBin(2);
                yticklabel = num2cell(bounds(1,1):diff(bounds(1,:)/5):bounds(1,2));
                xticklabel = num2cell(bounds(2,1):diff(bounds(2,:)/5):bounds(2,2));
                set(gca,'ytick',ytick,'yticklabel',yticklabel,'xtick',xtick,'xticklabel',xticklabel)
                ylabel(paramLabel{1}); xlabel(paramLabel{2})
                
                % add the estimated choice probability to the plot
                figure(trialfh(task,ct)); hold on;
                [foo foo foo leg] = legend;
                plot(poutEst{task,ct,2}.PP(:,1),'color',[0 .55 .55],'linewidth',2);
                leg{end+1} = 'est. p(choose blue)';
                legend(leg,'location','northeastoutside');
                
                % If simulating data, plot the correlation the recovered choice probabilities
                % and the actual choice probabilities
                if simulate;
                    psim = pout{task,ct}.PP(:,1);
                    hfig = figure;set(hfig,'position',[10 60 500 500],'Color','w'); box off;
                    [foo,foo,h(1)] = myScatter(psim,poutEst{task,ct,1}.PP(:,1),false,[0 0 1],'x');
                    [foo,foo,h(2)] = myScatter(psim,poutEst{task,ct,2}.PP(:,1),false,[1 0 0],'o');
                    h(3) = plot([0 1],[0 1],'k:','linewidth',2);
                    legend(h(1:2),{'Maximum Likelihood','Expected Value'},'location','best');legend boxoff
                    xlim([0 1]); ylim([0 1]);
                    xlabel('p(blue) - true'); ylabel ('p(blue) - estimated');
                    title('true vs. estimated choice probability');
                end
                
            end
            % - store the marginal likelihoods across subjects, sorted by task
            % - for each parameter, compute the peak (= maximum likelihood) and
            % the mean ('expected value', i.e. sum over all bins of (value(bin)
            % * likelihood(bin))
            
            fitted.alpha.pdf(task,ct,:) = marglik{1};
            fitted.beta.pdf(task,ct,:) = marglik{2};
            fitted.alpha.ml(task,ct) = ML(1);
            fitted.beta.ml(task,ct) = ML(2);
            fitted.alpha.ev(task,ct) = EV(1);
            fitted.beta.ev(task,ct) = EV(2);
        end % subject
        
        % PLOT posterior densities
        tick = 0:nBin/5:nBin;
        
        % plot posterior density of alpha across all subjects
        figure(fhParamPDF);
        subplot(2,2,task);hold on;
        plot(squeeze(fitted.alpha.pdf(task,:,:))');
        title(sprintf('%s volatility, alpha PDF',volatility{task}))
        xticklabel = num2cell(bounds(1,1):diff(bounds(1,:)/5):bounds(1,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{1})
        ylabel('p(data|m)');
        if simulate;	vline(simAlpha_bin,'k:'); end
        
        % plot posterior density of beta across all subjects
        subplot(2,2,task+2);hold on;
        plot(squeeze(fitted.beta.pdf(task,:,:))') ;
        title(sprintf('%s volatility, beta PDF',volatility{task}))
        xticklabel = num2cell(bounds(2,1):diff(bounds(2,:)/5):bounds(2,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{2})
        ylabel('p(data|m)');
        if simulate;	vline(simBeta_bin,'k:'); end
        
        
    end % task version
    
    if ~plotIndividual
        
        % plot barplots of the max and EV for alpha
        figure(fhParBar);
        
        subplot(2,2,1); barScatter(fitted.alpha.ml,[],volatility,true); box off;
        if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
        xlabel('volatility');
        ylabel('alpha (max. likelihood)');
        title('alpha - max likelihood');
        ylim([bounds(1,1) bounds(1,2)])
        
        subplot(2,2,2); barScatter(fitted.alpha.ev,[],volatility,true); box off;
        if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
        xlabel('volatility');
        ylabel('alpha (expected value)');
        title('alpha - expected value');
        ylim([bounds(1,1) bounds(1,2)])
        
        % plot barplots of the max and EV for beta
        subplot(2,2,3); barScatter(fitted.beta.ml,[],volatility,true); box off;
        if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
        xlabel('volatility');
        ylabel('beta (max. likelihood)');
        title('beta - max likelihood');
        ylim([bounds(2,1) bounds(2,2)])
        
        subplot(2,2,4); barScatter(fitted.beta.ev,[],volatility,true); box off;
        if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
        xlabel('volatility');
        ylabel('beta (expected value)');
        title('beta - expected value');
        ylim([bounds(2,1) bounds(2,2)])
        
        
        % add the trial-wise choice probabilities across subjects to the
        % probability plots
        if ~simulate
            figure(fhData)
            for task = 1:2
                nsub = length(sortedsubs{task});
                tmp = nan(nsub,nTrial);
                for isub = 1:nsub
                    tmp(isub,:)= poutEst{task,isub,2}.PP(:,1);
                end
                mu = mean(tmp,1);
                b = nan(nTrial,1,1);
                b(:,1,1) = std(tmp,[],1);
                subplot(2,1,task);hold on
                [foo foo foo leg] = legend;
                boundedline(1:nTrial,mu,b,'alpha','cmap',[0 .55 .55])
                leg{end+1} = 'estimated p(choice)';
                legend(leg,'location','northeastoutside');
                legend boxoff
            end
        end
    end
end

%==========================================================================
%% Section 4. Will be implemented in future versions
% Compute how good the model is: compute the likelihood, perhaps even plot
% trial-by-trial likelihood: what sort of trials can we explain well vs
% not? Can you think of other models and do model comparison?
%==========================================================================






