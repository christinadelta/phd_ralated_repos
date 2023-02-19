function fh = plotVFK_PVbin(volatility,learning_rate,predicted_state,probability,outc,pv)
% function fh = plotVFK_PVbin(volatility,learning_rate,predicted_state,probability,outc,pv)

% Date created: 18/02/2023

% plot stuff:
% volatility
% predected state
% for the 4 different combiations of lambda, v0 and omega values 

%% --------
% define plotting parameters, plotting font-size, font types, etc..
fnt         = 'Helvetica';
fn          = 'Arial';
abc         = 'ABCD';
fsl         = 14; 
fsA         = 18;
xsA         = -.15;
ysA         = 1.1;
yst         = 1.15; %0.15;
fst         = 18;
fsy         = 14;
% colours     = [1 .4 0;0 .4 1]; % orange and blue 
colours     = [0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980]; % green, orange
lcol        = [.8 .8 .8];

fpos0       = [0.2    0.0800    .85*1.0000    .65*0.8133];

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

% parameter names?
pnames      = {'\lambda','v_0','\omega'};
ylm         = [0 .59;-3 3];
nrows       = 2; % 2 rows of plots
ncols       = size(pv,1); % 4 columns 1 for each case of lambda, v0 and omega

% loop over combinations
for i = 1:size(pv,1)

    sub_plts = [i i+size(pv,1)]; % first value is for row 1, second value is for row 2

    if isrow(probability)
        dx = [0 diff(probability)~=0];
    else
        dx = [0; diff(probability)~=0];
    end

    ix = (find(dx));            % trial numbers where probabilities switch 
    tt = 1:length(probability); % trial indexes

    %% start with the plots 
    
    % plot volatility first
    fh(i,1) = subplot(nrows,ncols,sub_plts(1)); 
    plot(volatility(:,i),'color',colours(1,:),'linewidth',2); hold on;

    if ~any(isnan(ylm(1,:)))
        ylim(ylm(1,:)); 
    else
        ym = get(gca,'ylim');
        ym(2) = 1.01*ym(2);
        set(gca,'ylim',ym);
    end

    ym = get(gca,'ylim');

    % add gray vertical lines where probabilities change 
    for j=1:length(ix)
        plot([tt(ix(j)); tt(ix(j))],ym','color',lcol,'linewidth',2);
    end
    plot(volatility(:,i),'color',colours(1,:),'linewidth',2); hold on;
    set(gca,'fontname',fn);
    ylabel('Volatility','fontsize',fsy);
    set(gca,'ticklength', [0 0]);

    % plot predicted state
    fh(i,2) = subplot(nrows,ncols,sub_plts(2)); 
    plot(predicted_state(:,i),'color',colours(1,:),'linewidth',2); hold on;
    
    
    if ~any(isnan(ylm(2,:))), ylim(ylm(2,:)); end
    ym = get(gca,'ylim');

    % add gray vertical lines where probabilities switch
    for j=1:length(ix)
        plot([tt(ix(j)); tt(ix(j))],ym','color',lcol,'linewidth',2);
    end

    %vplot x and y labels
    plot(predicted_state(:,i),'color',colours(1,:),'linewidth',2); hold on;
    set(gca,'fontname',fn);
    ylabel('State predictions','fontsize',fsy); hold on;
    set(gca,'ticklength', [0 0]);
    xlabel('Trial','fontsize',fsy);

    % add case number 
    st = sprintf('Case %d',i);
    title(fh(i),st,'fontsize',fst,'fontname',fnt );


    % plot parameter values 
    st = sprintf('$%s=%0.1f$, $%s=%0.1f$, $%s=%0.1f$',pnames{1},pv(i,1),pnames{2},pv(i,2),pnames{3},pv(i,3));
    text(.5,1.12, st,'fontsize',16,'Unit','normalized','fontname',fnt,'parent',fh(i,2),...
         'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold','Interpreter','latex');

    % add plot number/label (ABCD)
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',fh(i,1));

end % end of combinations loop


end % end of function
