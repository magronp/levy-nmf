function [W,H,err] = cauchy_NMF(V,Wini,Hini,nNMF)

%init
W=Wini; H=Hini;
[F,T] = size(V);

% order of magnitude of W and H
W = W./ repmat(max(W,[],1),F,1);
H = H./ repmat(max(H,[],2),1,T);
H = max(V(:))*H;
WH = W*H;

%cost
err = zeros(1,nNMF+1);
err(1) = cost_cauchy(WH,V);

for it=1:nNMF
    
    a = 3/4 * (WH ./ (WH.^2+V.^2+eps)) * H';
    b = (WH+eps).^(-1) * H';
    
    W = W .* (b ./ (a+sqrt(a.^2+2*b.*a)+eps));

    WH = W*H;
    a = 3/4 * W'* (WH ./ (WH.^2+V.^2+eps)) ;
    b = W' * (WH+eps).^(-1) ;
    
    H = H .* (b ./ (a+sqrt(a.^2+2*b.*a)+eps));
 
    WH = W*H;
    
    % Normalization L1
    sumW2 = sum(W);
    W = W * diag(1./sumW2);
    H = diag(sumW2) * H;
    
    % cost
    err(it+1) = cost_cauchy(WH,V);
    
end

end


function d = cost_cauchy(X,Y)

X=X+eps; Y=Y+eps;
D=3/2*log(X.^2+Y.^2)-log(X);
d = sum(D(:));

end