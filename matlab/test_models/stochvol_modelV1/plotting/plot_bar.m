function h = plot_bar(nr,nc,sub_plts,mx,ex,labels,mstr,abc,col,yls)

if ~iscell(sub_plts)
    sub_plts = num2cell(sub_plts);
end

if nargin<8, abc = ''; end
if nargin<9, col = []; end
if nargin<10, yls = cell(1,length(mx)); end

fs = def_actions('fs');
fn = def_actions('fn');
fsy = def_actions('fsy');
alf = def_actions('alf');
fsA = def_actions('fsA');
xsA = def_actions('xsA');
ysA = def_actions('ysA');

if isempty(abc)
    abc = [];    
end

if isempty(col)
    col = def_actions('col');
end

if isempty(ex)
    ex = cell(1,length(mx));
    for i=1:length(mx)
        ex{i} = zeros(size(mx{i}));
    end    
end

bw = .025*numel(mx{1});

h = nan(1,length(mx));
for i=1:length(mx)
    h(i) = subplot(nr,nc,sub_plts{i});    
    errorbarKxN(mx{i}',ex{i}',labels,col,bw);    
    set(gca,'fontsize',fs);
    alpha(alf);
    ylabel(mstr{i},'fontsize',fsy);
    ax = ancestor(gca, 'axes');
    xaxes = get(ax,'XAxis');
    set(xaxes,'fontsize',fsy);    
    
    if ~isempty(abc)
        text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn);
    end
    if ~isempty(yls{i})
        ylim(yls{i});
    end
end

end
