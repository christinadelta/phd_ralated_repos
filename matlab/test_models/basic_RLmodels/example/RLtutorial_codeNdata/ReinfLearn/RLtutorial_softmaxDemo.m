function RLtutorial_softmaxDemo
% FUNCTION RLtutorial_softmaxDemo
% 
% short demo to show the effect of the value of beta in the softmax
% function
%
% The code computes choice probabilities based on the value difference
% between A and B, where vA = 1-vB, and vA ranges from 0-1, using the
% softmax function:
% 
%    p(A) = exp(beta*vA)/(exp(beta*vA)_exp(beta*vB))
% 
% and then generates a choice based on these probabilities, generating a
% random number r in the [0 1] range, where  if r<=p(A), then choice  = A,
% and B otherwise. 
% 
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------
 
betaArray = [3 9 10];

if any(abs(betaArray)>700)
    error('Please pick beta values between -700 and +700, as matlab cannot deal with such small/big numbers!');
end


h = figure; set(h,'position',[10 60 400 400 ],'Color','w'); 
hold on; box off; 

va = 0:0.02:1;
vb = 1-va; % value of B is [1 - value of A]
pa = nan(length(betaArray),length(va));
data = nan(length(betaArray),length(va));
ct = 0;
x = .02
for beta = betaArray
    ct = ct+1;
    pa(ct,:) = exp(beta*va)./(exp(beta*va)+exp(beta*vb)); % probability of choosing A
    leg{ct} = sprintf('beta = %d',beta);
    tmp = double(rand(1,length(va))< pa(ct,:))*(1+x*ct); % generate choices based on this data
    tmp(tmp==0)=tmp(tmp==0)-x*ct;
    data(ct,:) = tmp;
end

plot(va-vb,pa);
legend(leg,'location','best'); legend boxoff;
plot(va-vb,data,'*');
xlabel('Value(A) - Value (B)');
ylabel('p(choice=A)');
title(sprintf('softmax demo'));
ylim([-x*ct 1+x*ct]);

