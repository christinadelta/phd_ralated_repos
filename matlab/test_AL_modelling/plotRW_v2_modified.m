function h = plotRW_v1_modified(Qvals, Ps, a, x)
% Plot Q values and choice probabilities for each condition separately, along with the reward rate (x)

% Define plotting parameters, etc.
% (Initialization and figure setup code remains the same)

% Assuming x is structured similarly to Qvals and Ps, adjust as needed
nStc = 3; % Number of conditions
nVol = 2; % Number of volatilities per condition
h = zeros(nStc * nVol, 1); % Preallocate for subplot handles

% Iterate over conditions and volatilities
for ss = 1:nStc
    for vv = 1:nVol
        % Calculate subplot index
        subplotIndex = (ss-1) * nVol + vv;
        
        % Select subplot and store its handle
        h(subplotIndex) = subplot(nStc, nVol, subplotIndex);
        
        % Extract data for current condition and volatility
        currentQvals = Qvals{ss, vv};
        currentPs = Ps{ss, vv};
        currentA = a{ss, vv};
        currentX = x{ss, vv}; % Extract reward rate for current condition and volatility
        
        % Plotting commands
        plot(currentQvals, 'LineWidth', 2);
        hold on;
        plot(currentPs, '--', 'LineWidth', 2);
        plot(currentX, '-.', 'LineWidth', 2); % Plot reward rate
        
        % Adjust plot (titles, labels, legends)
        title(sprintf('Condition %d, Volatility %d', ss, vv));
        xlabel('Trials');
        ylabel('Value/Probability/Reward Rate');
        legend({'Q values', 'Choice probabilities', 'Reward Rate'}, 'Location', 'best');
        
        hold off;
    end
end

% Adjust figure properties as necessary
set(gcf, 'Name', 'Q values, Choice Probabilities, and Reward Rate by Condition', 'NumberTitle', 'off');
end
