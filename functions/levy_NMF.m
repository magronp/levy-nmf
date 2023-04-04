%  Lévy NMF algorithm under Multiplicative Update Rules
%
%  Refeferences:
%  "Lévy NMF for robust nonnegative source separation", P. Magron, R.Badeau and A. Liutkus, Proc. IEEE WASPAA, Oct. 2017
%  "Lévy NMF : un modèle robuste de séparation de sources non-négatives", P. Magron, R.Badeau and A. Liutkus, Proc. GRETSI, Sept. 2017
%
%  Code by Paul Magron, march 2017.
% 
% Inputs :
%     V : nonnegative data matrix F*T
%     Wini : initial dictionary matrix F*K
%     Hini : initial activation matrix K*T
%     nNMF : number of iterations
% 
% Outputs :
%     W : estimated dictionary
%     H : estimated activation matrix
%     err : cost function

function [W,H,err] = levy_NMF(V,Wini,Hini,nNMF)

%init
W=Wini; H=Hini;
[F,T] = size(V);

% order of magnitude of W and H
W = W./ repmat(max(W,[],1),F,1);
H = H./ repmat(max(H,[],2),1,T);
H = sqrt(max(V(:)))*H;
Vhat = W*H;

%cost
err = zeros(1,nNMF+1);
err(1) = ISdiv(Vhat.^2,V);

for it=1:nNMF
    
    % up W
    W = W .* ( ( ((Vhat+eps).^(-1))*H' )  ./  ((Vhat.*(V+eps).^(-1))*H' +eps) ).^0.5;
    
    Vhat = W*H;
    
    % up H
    H = H .* ( ( W'*((Vhat+eps).^(-1)) )  ./  (W'*(Vhat.*(V+eps).^(-1)) +eps) ).^0.5;
 
    % Normalization L1
    sumW = sum(W);
    W = W * diag(1./sumW);
    H = diag(sumW) * H;
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
