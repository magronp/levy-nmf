clc; clear all; close all;
Fs = 44100; Nw = 4096; Nfft = 4096; hop = Nw/4;
nNMF=50;
Nbins = 1000;
nmax = 5;

% Signal generation
source_type='SYNTH_NO_OL';
gen_sources_time; gen_sources_TF;

% Initialization
V = abs(X);
Wini = rand(F,K); Hini = rand(K,T);

% Levy NMF
[Wl,Hl] = levy_NMF(V,Wini,Hini,nNMF);
Se_levy = wiener(V,Wl,Hl);

% KLNMF : Poisson model
[Wp,Hp] = NMF(V,Wini,Hini,nNMF,1,0);
Se_poisson = wiener(V,Wp,Hp);

% Histograms
histo_levy = zeros(K,Nbins); edges_levy = zeros(K,Nbins+1);
for k=1:K
    auxSe = Se_levy(:,:,k) ./ ((Wl(:,k)*Hl(k,:)).^2+eps); auxSe = auxSe(:);
    aux = auxSe(auxSe<nmax);
    [histo_levy(k,:),edges_levy(k,:)]= histcounts(aux,Nbins);
    histo_levy(k,:)=histo_levy(k,:)./max(histo_levy(k,:));
end
bins_levy = edges_levy(:,1:end-1)+nmax/Nbins;

scale_levy = zeros(1,K);
[~,ind] = max(histo_levy,[],2);
for k=1:K
    scale_levy(k) = 3*bins_levy(k,ind(k));
end

histo_poisson = zeros(K,Nbins); edges_poisson = zeros(K,Nbins+1);
for k=1:K
    auxSe = Se_poisson(:,:,k)./ ((Wl(:,k)*Hl(k,:)).^2+eps); auxSe = auxSe(:);
    aux = auxSe(auxSe<nmax);
    [histo_poisson(k,:),edges_poisson(k,:)]= histcounts(aux,Nbins);
    histo_poisson(k,:)=histo_poisson(k,:)./max(histo_poisson(k,:));
end
bins_poisson = edges_poisson(:,1:end-1)+nmax/Nbins;
scale_poisson = zeros(1,K);
[~,ind] = max(histo_poisson,[],2);
for k=1:K
    scale_poisson(k) = bins_poisson(k,ind(k));
end

% Plot
k=1;
lnn=levy_density(bins_levy(k,:),scale_levy(k));
pnn=poisson_density(bins_poisson(k,:),scale_poisson(k));

figure;
subplot(1,2,1);
bar(bins_levy(k,:),histo_levy(k,:).*max(lnn)); hold on;
plot(bins_levy(k,:),lnn,'r');

subplot(1,2,2);
bar(bins_poisson(k,:),histo_poisson(k,:).*max(pnn)); hold on;
plot(bins_poisson(k,:),pnn,'r');


% Synthesis
xel = zeros(size(sm)); xep = zeros(size(sm));
for k=1:K
    xel(k,:) = iSTFT(Se_levy(:,:,k).*exp(1i*angle(X)),Nfft,w,hop);
    xep(k,:) = iSTFT(Se_poisson(:,:,k).*exp(1i*angle(X)),Nfft,w,hop);
end


