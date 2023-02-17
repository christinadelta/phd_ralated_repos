function [o, x] = timeseries_bin


simcat = 'basic';
p0 = .8;
p  = [p0*ones(1,2) repmat([1-p0 p0],1,2) (1-p0)*ones(1,4)];
n  = 20;

fname = 'bin.mat';

pipedir = getdefaults('pipedir');
fdir = fullfile(pipedir,simcat); makedir(fdir);
fname = fullfile(fdir,fname);

if ~exist(fname,'file')
    nb = length(p);    
    N  = nb*n;
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