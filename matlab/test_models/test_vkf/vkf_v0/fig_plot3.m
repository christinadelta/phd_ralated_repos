function fig_plot3(x,y,val1,vol1,lr1,val2,vol2,lr2)
% plot vkf and hgf results (figure 6)

nr = 3;
nc = 2;
fpos0 = [0.2    0.0800    .55*1.0000    .65*0.8133];


fn = getdefaults('fn');
fnt = getdefaults('fnt');
fst = getdefaults('fst');
% fsl = getdefaults('fsl');
fsA = getdefaults('fsA');
xsA = getdefaults('xsA');
ysA = getdefaults('ysA');
abc = getdefaults('abc');
yst = getdefaults('yst');

%---------
figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

yl = [0 1.5; .2 1.2;-2.5 2.5];


sub_plts = [1 3 5];
hl = sim_C_plot(nr,nc,sub_plts,x,y,vol1,lr1,val1,yl);

text(.5,yst,'VKF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hl(1),'HorizontalAlignment','Center','fontweight','bold');


for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hl(i));
end
%---------
sub_plts = [2 4 6];

hr = sim_C_plot(nr,nc,sub_plts,x,y,vol2,lr2,val2,yl);
text(.5,yst,'HGF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hr(1),'HorizontalAlignment','Center','fontweight','bold');

for i=1:2
    text(xsA,ysA,abc(i+2),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hr(i));
end
end