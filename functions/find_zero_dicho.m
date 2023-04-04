function [xm,xp] = find_zero_dicho(a,Niter)

[N,M] = size(a);

cost_diff=@(eta) exp(2.*eta.*log(a))-1-2.*eta.*a.*log(a);

% Initialize with 1/2 and Newton method (to ensure to be > than the zero)
xm0=1/2*ones(N,M);
eta0=ones(N,M);
xp0 = eta0+(1+2.*a.*log(a).*eta0-a.^(2.*eta0))./(2.*a.*log(a).*(a.^(2.*eta0-1)-1));

xm = xm0; xp=xp0;
c=0;

% Dichotomy loop
while (c<Niter)
    xmed = (xm+xp)/2;
    costmed = cost_diff(xmed);
    
    xm(costmed<=0)=xmed(costmed<=0);
    xp(costmed>0)=xmed(costmed>0);
    
    c=c+1;
end

end
