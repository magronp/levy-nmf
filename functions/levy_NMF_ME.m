function [W,H,err] = levy_NMF_ME(V,Wini,Hini,nNMF,mask,algo,meth)

%init
W=Wini; H=Hini;
[F,T] = size(V);

if nargin<7
    meth = 'DL';
end

if nargin<6
    algo = 'MM';
end

if nargin<5
    mask = ones(F,T);
end

% order of magnitude of W and H
% W = W./ repmat(max(W,[],1),F,1);
% H = H./ repmat(max(H,[],2),1,T);
% H = sqrt(max(V(:)))*H;
Vhat = W*H;

%cost
err = zeros(1,nNMF+1);
err(1) = ISdiv(Vhat.^2,V);

for it=1:nNMF
    
    %Update W
    aW = ( ( mask.*((Vhat+eps).^(-1))*H' )  ./  ((mask.*Vhat.*(V+eps).^(-1))*H' +eps) );
    
    switch algo
        case 'H'
            etaW=1;
        case 'MM'
            etaW=1/2;
        case 'ME'
            etaW=compute_eta_opt(aW,meth);
    end
    
    W = W .* (aW+eps).^(etaW+eps);
    
    Vhat = W*H;
    
    %Update H
    aH =( ( W'*(mask.*(Vhat+eps).^(-1)) )  ./  (W'*(mask.*Vhat.*(V+eps).^(-1)) +eps) );
    switch algo
        case 'H'
            etaH=1;
        case 'MM'
            etaH=1/2;
        case 'ME'
            etaH=compute_eta_opt(aH,meth);
    end
    
    H = H .* (aH+eps).^(etaH+eps);
 
    % Normalization L1
    sumW2 = sum(W);
    W = W * diag(1./sumW2);
    H = diag(sumW2) * H;
    Vhat = W*H;
    
    % cost
    err(it+1) = ISdiv(Vhat.^2,V);
    
end

end


function d = ISdiv(X,Y)

X=X+eps; Y=Y+eps;
E = X./Y - log(X./Y) -1;
d = sum(E(:));

end

function eta = compute_eta_opt(a,meth)
             
switch meth
    
    case 'line'
        eta= 1/2 + (a.*log(a)-a+1)./(2*(a.^2-a-a.*log(a))+eps); % corde 1/2->1 
        
    case 'DL'
        eta = 1/2+(a.*log(a)-a+1)./(a-1).^2; %DL order 1
        
    case 'dicho'
        eta = find_zero_dicho(a,10);
        
    case 'newton'
        eta=ones(size(a));
        for iter=1:10
            eta = eta+(1+2.*a.*log(a).*eta-a.^(2.*eta))./(2.*a.*log(a).*(a.^(2.*eta-1)-1));
        end
        %eta = 1+(1+2.*a.*log(a).*1-a.^(2.*1)+eps)./(2.*a.*log(a).*(a.^(2.*1-1)-1)+eps);

end

eta(a<=1) = 1/2;

end