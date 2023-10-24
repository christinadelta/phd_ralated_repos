function [vol,stc,lr,val,unc] = model_parfilter(o,specs,islesioned)
if nargin<3, islesioned = 'healthy'; end

np = specs.nparticles;
% np = 100;
% % 
% % specs.lambda_v  = params(1);
% % specs.lambda_s  = params(2);
% % specs.v0        = params(3);
% % specs.s0        = params(4);
% specs.x0_unc    = 1;


x0_unc = 1;
if isfield(specs,'x0_unc')
    x0_unc = specs.x0_unc;
end    

switch lower(islesioned(1:3))
    case 'hea'
        lambda_v = specs.lambda_v;
        lambda_s = specs.lambda_s;
        v = specs.v0;
        s = specs.s0;
        
%         lambda_v    = params(1);
%         lambda_s    = params(2);
%         v           = params(3);
%         s           = params(4);

        [val,vol,stc,lr,unc] = parfilter_core(o,x0_unc,lambda_v,lambda_s,v,s,np);   

    case 'vol'
        lambda_s = specs.lambda_s;
        v_lesioned = specs.v0_lesioned;        
        s = specs.s0;
        [val,vol,stc,lr,unc] = parfilter_core_vol_lesioned(o,x0_unc,lambda_s,v_lesioned,s,np);        
    case 'sto'
        lambda_v = specs.lambda_v;
        s_lesioned = specs.s0_lesioned;
        v = specs.v0;
        [val,vol,stc,lr,unc] = parfilter_core_stc_lesioned(o,x0_unc,lambda_v,v,s_lesioned,np);        
    otherwise
        error('bad 3rd input: %d',islesioned);
end

