function Xw = wiener(X,W,H)

[F,T] = size(X);
[~,K] = size(W);
Xw = zeros(F,T,K);
for k=1:K
  Xw(:,:,k) = W(:,k) * H(k,:);
end

end