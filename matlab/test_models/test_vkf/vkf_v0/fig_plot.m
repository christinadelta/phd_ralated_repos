function fig_plot(x,v,lr,m,y)

fpos0 = [0.2    0.0800    .55*1.0000    .7*0.8133];


fn = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
fsl = getdefaults('fsl');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');

yst = getdefaults('yst');

%---------

figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

%------------
nr = 3;
nc = 2;
% sub_plts = [1 3 5];
% 
% yl = [0 .4;.5 .85;-1.9 1.9];

% j = 1;
% h = sim_A_plot(nr,nc,sub_plts,x(:,j),y(:,j),v(:,j),lr(:,j),m(:,j),yl);
% 
% for i=1:3
%     text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
% %     set(h(i),'ylim',yl(i,:));
% end
% 
% i =1;
% text(.5,yst,'Linear','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');

%------------
% sub_plts = [2 4 6];
sub_plts = [1 2 3];

yl = [.08 .16;.35 .5;-.2 1.2];

% j = 2;
j = 1;
h = sim_A_plot(nr,nc,sub_plts,x(:,j),y(:,j),v(:,j),lr(:,j),m(:,j),yl);

for i=1:3
    text(xsA,ysA,abc(i+3),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end
hlg = legend(h(3),{'Predicted','True'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);

i =1;
text(.5,yst,'Binary','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');

end