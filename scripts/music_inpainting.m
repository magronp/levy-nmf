%%% Restoration of music spectrograms corrupted with impulsive noise

clear all; clc; close all;
pkg load signal statistics;

% Random seed for reproducibility
rand("state", 1);

% Data parameters
Fs = 8000;
Ndata = 6;
dataset_path = 'data/guitar/';

% STFT parameters
Nfft = 1000; Nw = 1000; hop = Nw/4; wtype = 'hann';

% NMF parameters
Nnmf = 200;
K = 30;

Ndata=6;
score = zeros(Ndata,6);
KLdiv = zeros(Ndata,6);

for it=1:Ndata
  
    clc; fprintf(' data %d / %d \n',it,Ndata);
    
    [x,X,F,T,ts,freq] = get_data_guitar(dataset_path,it,Fs,Nfft,hop,Nw,wtype);
    V = abs(X);
    
    % Initial NMF matrices
    Wini = rand(F,K); Hini = rand(K,T);

    % Corruption of the spectrograms
    corr = rand(size(V)); corr(corr<0.9) = 0;
    V_corr = max(V,max(V(:))*corr);
    
    %With or without localization of the noise
    %delta = 1-corr; delta(delta<1)=0;
    delta = ones(F,T);
    
    % Reestimation
    [Wis_rest,His_rest] = NMF(V_corr.^2,Wini,Hini,Nnmf,0,0,0,delta);
    Vis_rest = sqrt(Wis_rest*His_rest);

    [Wl_rest,Hl_rest,err] = levy_NMF(V_corr,Wini,Hini,Nnmf);
    Vl_rest = (Wl_rest*Hl_rest).^2;

    [Wkl_rest,Hkl_rest] = NMF(V_corr,Wini,Hini,Nnmf,1,0,0,delta);
    Vkl_rest = Wkl_rest*Hkl_rest;
    
    [Wc_rest,Hc_rest] = cauchy_NMF(V_corr,Wini,Hini,Nnmf);
    Vc_rest = Wc_rest*Hc_rest;
    
    aux = inexact_alm_rpca(V_corr);
    Vr_rest = abs(aux);
    
    % Synthesis
    Xcorr = V_corr.* exp(1i*angle(X)); xcorr = iSTFT(Xcorr, Nfft, hop, Nw, wtype);
    Xl = Vl_rest .* exp(1i*angle(X)); xl =iSTFT(Xl, Nfft, hop, Nw, wtype);
    Xis = Vis_rest .* exp(1i*angle(X));  xis =  iSTFT(Xis, Nfft, hop, Nw, wtype);
    Xc = Vc_rest .* exp(1i*angle(X));  xc =  iSTFT(Xc, Nfft, hop, Nw, wtype);
    Xr = Vr_rest .* exp(1i*angle(X));  xr = iSTFT(Xr, Nfft, hop, Nw, wtype);
    Xkl = Vkl_rest .* exp(1i*angle(X));  xkl = iSTFT(Xkl, Nfft, hop, Nw, wtype);

    % Score
    score(it,1) = GetSDR(xcorr,x);
    score(it,2) = GetSDR(xis,x);
    score(it,3) = GetSDR(xkl,x);
    score(it,4) = GetSDR(xc,x);
    score(it,5) = GetSDR(xr,x);
    score(it,6) = GetSDR(xl,x);
    
    % KL divergence between spectro
    KLdiv(it,1) = beta_div(V_corr,V,1);
    KLdiv(it,2) = beta_div(Vis_rest,V,1);
    KLdiv(it,3) = beta_div(Vkl_rest,V,1);
    KLdiv(it,4) = beta_div(Vc_rest,V,1);
    KLdiv(it,5) = beta_div(Vr_rest,V,1);
    KLdiv(it,6) = beta_div(Vl_rest,V,1);

end

% Plot spectrograms
limtsplot = [-4 max(log10(V_corr(:)))];

figure;
subplot(2,3,1); imagesc(ts,freq,log10(V),limtsplot); axis xy; ylabel('Fréquence (Hz)','fontsize',16); title('Original','fontsize',16);
subplot(2,3,2); imagesc(ts,freq,log10(V_corr),limtsplot); axis xy;  title('Corrompu','fontsize',16);
subplot(2,3,3); imagesc(ts,freq,log10(Vkl_rest),limtsplot); axis xy; title('ISNMF','fontsize',16);
subplot(2,3,4); imagesc(ts,freq,log10(Vc_rest),limtsplot); axis xy; xlabel('Temps (s)','fontsize',16); ylabel('Fréquence (Hz)','fontsize',16); title('Cauchy NMF','fontsize',16);
subplot(2,3,5); imagesc(ts,freq,log10(Vr_rest),limtsplot); axis xy; xlabel('Temps (s)','fontsize',16); title('RPCA','fontsize',16);
subplot(2,3,6); imagesc(ts,freq,log10(Vl_rest),limtsplot); axis xy; xlabel('Temps (s)','fontsize',16); title('Lévy NMF','fontsize',16);
hc=colormap(gray); hc=hc(end:-1:1,:); colormap(hc);

% Plot results (KL div and SDR)
algos = {'Corrupted','ISNMF','KLKNMF','Cauchy NMF','RPCA','Lévy NMF'};
Nalgos = length(algos);

figure; boxplot(log10(KLdiv));
set(gca,'xtick', 1:Nalgos, 'xticklabel', algos,'fontsize',14);
h=ylabel('$\log (KL)$'); set(h,'Fontsize',16,'interpreter','latex');

figure; boxplot(score);
set(gca,'xtick', 1:Nalgos, 'xticklabel', algos,'fontsize',14);
h=ylabel('SDR (dB)'); set(h,'Fontsize',16);
