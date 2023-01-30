function figh = plot_manysubs(allQs, allPs, choices, fprob, params, condition)

% created on 30/01/2023
% plot averaged (many simulated subs)

% -----------------------------

% make figure handle
figh        = figure('Name','Data'); box off; hold on;
set(figh,'position',[10 60 900 650],'paperunits','centimeters','Color','w');
figure(figh)

if condition == 1

    % average q values, choices and choice probabilities over subjects 
    avQs_v      = mean(squeeze(allQs(:,:,1)),1);
    avQs_h      = mean(squeeze(allQs(:,:,2)),1);
    avPs_v      = mean(allPs,1);
    avchoices   = nanmean(choices,1);

    % make plot
    plot(fprob, 'k:', 'LineWidth', 2) % plot the fedback probabilities
    plot(avchoices, 'k-', 'LineWidth',1)
    plot(avPs_v, 'b-', 'LineWidth', 2)
    plot(avQs_v, ':', 'Color',[.4 0.4 1], 'LineWidth',2)
    plot(avQs_h, ':', 'Color',[1 .35 0.1], 'LineWidth',2)
    legend({'p(reward|verical)','mean choice (vertical)', 'p(choose vertical)','value(vertical)','value(horizontal)',...
    },'location','northeastoutside');

    legend boxoff
    ylabel('probability');
    ylim([-0.1 1.1]);
    xlabel('trial');
    title(sprintf('mean data, alpha = %0.02f, beta = %02.01f',...
        params(1),params(2)));

else % if volatility is added 

    
        



end % end of condition statement 


end % end of function