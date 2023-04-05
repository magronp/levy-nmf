
function D = beta_div(X,Y,beta,mask)

X=X+eps; Y=Y+eps;

if nargin<4
    mask = ones(size(X));
end

switch beta
    case 0
        E = X./Y - log(X./Y) -1;
    case 1
        E = X .* log(X./Y)+Y-X;
    otherwise
        E =1/(beta*(beta-1)) * (X.^beta + (beta-1)*Y.^beta - beta*X.*Y.^(beta-1));
end

Emask = E .* mask;
D = sum(Emask(:));

end