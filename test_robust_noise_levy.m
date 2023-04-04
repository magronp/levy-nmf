clear all; clc; close all;
Fs = 8000; Nfft = 1000; Nw = 1000; hop = Nw/4;
Nnmf=200;
K=30;
Ndata=1;
score = zeros(Ndata,6);
KLdiv = zeros(Ndata,6);

for it=1:Ndata
    clc; fprintf(' data %d / %d \n',it,Ndata);
    num_piece = it;
    [x,X,w,F,T,ts,freq] = get_data_guitar_piece(Fs,Nfft,hop,num_piece);
    
    V = abs(X);
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
    Xcorr = V_corr.* exp(1i*angle(X)); xcorr = iSTFT(Xcorr,Nfft,hop);
    Xl = Vl_rest .* exp(1i*angle(X)); xl = iSTFT(Xl,Nfft,hop);
    Xis = Vis_rest .* exp(1i*angle(X));  xis = iSTFT(Xis,Nfft,hop);
    Xc = Vc_rest .* exp(1i*angle(X));  xc = iSTFT(Xc,Nfft,hop);
    Xr = Vr_rest .* exp(1i*angle(X));  xr = iSTFT(Xr,Nfft,hop);
    Xkl = Vkl_rest .* exp(1i*angle(X));  xkl = iSTFT(Xkl,Nfft,hop);

    % Score
    score(it,1) = bss_eval_sources(xcorr,x);
    score(it,2) = bss_eval_sources(xis,x);
    score(it,3) = bss_eval_sources(xkl,x);
    score(it,4) = bss_eval_sources(xc,x);
    score(it,5) = bss_eval_sources(xr,x);
    score(it,6) = bss_eval_sources(xl,x);
    
    % KL divergence between spectro
    KLdiv(it,1) = beta_div(V_corr,V,1);
    KLdiv(it,2) = beta_div(Vis_rest,V,1);
    KLdiv(it,3) = beta_div(Vkl_rest,V,1);
    KLdiv(it,4) = beta_div(Vc_rest,V,1);
    KLdiv(it,5) = beta_div(Vr_rest,V,1);
    KLdiv(it,6) = beta_div(Vl_rest,V,1);

end


% Scaling for plots
% min_plot = 10^(-1);
% V = max(V,min_plot); V_corr = max(V_corr,min_plot);
% Vis_rest = max(Vis_rest,min_plot);
% Vkl_rest = max(Vkl_rest,min_plot);
% Vc_rest = max(Vc_rest,min_plot);
% Vr_rest = max(Vr_rest,min_plot);
% Vl_rest = max(Vl_rest,min_plot);

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
%colormap(jet);

% Display KL divergence
algos = {'Corrompu','ISNMF','KLKNMF','Cauchy NMF','RPCA','Lévy NMF'};
figure; boxplot(log10(KLdiv),algos);
h=ylabel('$\log (KL)$'); set(h,'Fontsize',16,'interpreter','latex');

% Display DSR
figure; boxplot(score,algos);
h=ylabel('SDR (dB)'); set(h,'Fontsize',16);

