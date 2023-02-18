function fig_plot4(x,v,m,tx)

fpos0 = [0.2    0.0800    .85*1.0000    .65*0.8133];

fn  = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');
%---------

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

%------------
pnames = {'\lambda','v_0','\omega'};

yl = [0 .59;-3 3];

nr = 2;
nc = size(m,2);

for j = 1:size(m,2)
    sub_plts = [j j+size(m,2)];
    
    h(j,:) = sim_plot(nr,nc,sub_plts,x,v(:,j),m(:,j),yl); %#ok<AGROW>
    
    st = sprintf('Scenario %d',j);
    title(h(j),st,'fontsize',fst,'fontname',fnt );
    
    st = sprintf('$%s=%0.1f$, $%s=%0.1f$, $%s=%0.1f$',pnames{1},tx(j,1),pnames{2},tx(j,2),pnames{3},tx(j,3));
    text(.5,1.12, st,'fontsize',16,'Unit','normalized','fontname',fnt,'parent',h(j,2),...
         'HorizontalAlignment','center','VerticalAlignment','middle','fontweight','bold','Interpreter','latex');
    
    for i=1:2
    end     
    text(xsA,ysA,abc(j),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(j,1));
    
end

end
