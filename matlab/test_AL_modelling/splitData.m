function blockdata = splitData(stc,vol,x,alldata,stdvals)

%  PAL_PROJECT
% the function split all the data required for fitting into
% blocks/conditions, for the kalman model!

%% 

% how many subjects do we have?
subs = size(alldata,1);

% for each subject split all data (outcomes, actions,rts and x) across
% conditions/blocks. Also apply GNA to the outcomes 
for sub = 1:subs

    a   = alldata{sub,1}.actions;
    o   = alldata{sub,1}.outcome;
    rt  = alldata{sub,1}.rts;

    
    o(find(o==2))   = 0; % change outcomes into 1 and 0
    oR              = 1 - o; % we need to compute values for both outcomes 

    for ss = 1:3 % loop over stc levels

        tmpx    = x(stc(:,ss),:);
        tmpo    = o(stc(:,ss),:);
        tmpa    = a(stc(:,ss),:);
        tmprt   = rt(stc(:,ss),:);
        tmpoR   = oR(stc(:,ss),:);

        % add GNA to outcomes (for modelling using the kalman filter
        std_o   = tmpo + stdvals(ss) * randn(size(tmpo));
        std_oR  = tmpoR + stdvals(ss) * randn(size(tmpoR));

        for t = 1:2 % loop over vol phases

            blockdata{1,sub}{ss,t}.x     = tmpx(vol(:,t),:);
            blockdata{1,sub}{ss,t}.o     = std_o(vol(:,t),:);
            blockdata{1,sub}{ss,t}.oR    = std_oR(vol(:,t),:);
            blockdata{1,sub}{ss,t}.a     = tmpa(vol(:,t),:);
            blockdata{1,sub}{ss,t}.rt    = tmprt(vol(:,t),:);

        end % end of vol loop
       
    end % end of ss loop
end % end of subs loo

end % end of function