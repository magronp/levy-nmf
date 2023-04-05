function [W,H] = sort_distcorr(Wref,W,H)

K=size(Wref,2);
co = zeros(K,K);

for k1=1:K
    for k2=1:K
        co(k1,k2) = distcorr(Wref(:,k1),W(:,k2));
    end
end

permut = zeros(K,K);
aux = co;

for k=1:K
    [~,ind] = max(aux(:));
    [ir,ic] = ind2sub(size(aux),ind);
    permut(ir,ic)=1;
    aux(ir,:) = 0; aux(:,ic)=0;
end

W = W * permut';
H = permut * H;


end