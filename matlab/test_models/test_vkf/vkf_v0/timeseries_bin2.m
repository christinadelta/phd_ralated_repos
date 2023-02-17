function [o, x] = timeseries_bin2
simcat = 'basic';
n = 10;
b = .8;
fname = sprintf('bin_n%d.mat',n);

p = [b*ones(1,2) repmat([1-b b],1,2) (1-b)*ones(1,4)];
p = [1-p p 1-p p];

pipedir = getdefaults('pipedir');
fdir = fullfile(pipedir,simcat); makedir(fdir);
fname = fullfile(fdir,fname);
if ~exist(fname,'file')        
    nb = length(p);  
    N  =  nb*n;    
    x  = nan(1,N);    
    o  = zeros(1,N);

    t0 = 0;
    for i=1:nb
        ii = t0 + (1:n);
        x(ii) = p(i);
        ni = randperm(n);
        ni = ni(1: round(p(i)*n));
        o(ii(ni))  = 1;
        t0 = t0 + n;
    end
    
    save(fname,'o','x');
end
data = load(fname);
o = data.o';
x = data.x';

end
