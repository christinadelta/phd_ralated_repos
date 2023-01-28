function trl_plts = trlplots(avtrl_choices, alphas, betas)

% initiate figure handle
trl_plts = figure; box off; hold on;

% plot alpha = 0.2
subplot(2,2,1)
alpha1 = plot(avtrl_choices{1,1}');
xlabel('trials')
ylabel('choice p(vertical)')
title([' learning  curve \alpha = ', num2str(alphas(1))])

% plot alpha = 0.4
subplot(2,2,2)
alpha3 = plot(avtrl_choices{1,2}');
xlabel('trials')
ylabel('choice p(vertical)')
title([' learning  curve \alpha = ', num2str(alphas(2))])

% plot alpha = 0.7
subplot(2,2,3)
alpha3 = plot(avtrl_choices{1,3}');
xlabel('trials')
ylabel('choice p(vertical)')
title([' learning  curve \alpha = ', num2str(alphas(3))])

% plot alpha = 1
subplot(2,2,4)
alpha4 = plot(avtrl_choices{1,4}');
xlabel('trials')
ylabel('choice p(vertical)')
title([' learning  curve \alpha = ', num2str(alphas(4))])

% make legend 
for i = 1:length(betas)

    leg{i} = ['\beta = ' num2str(betas(i))];

end

leg1 = legend(alpha1(end:-1:1), {leg{end:-1:1}});

set(leg1, 'fontsize', 12)
%set(leg1, 'position', [0.6267    0.6453    0.1500   0.2000]);


end
