clear all; clc; close all;
Fs = 44100; Nfft = 4096; Nw = 4096; hop = Nw/4;
source_type='DSD100';
Nnmf=200;
direc = 'levy NMF/sounds/';

Ndata=1;
score = zeros(Ndata,5);
score2 = zeros(Ndata,5);

for it=1:Ndata
    clc; fprintf(' iter %d / %d \n',it,Ndata);
    num_piece = it;
    gen_sources_time; gen_sources_TF_inst; K=30;

    smusic = sum(sm(1:3,:),1);
    svoice = sm(4,:);
    
    Vmusic = abs(sum(Sm(:,:,1:3),3));
    Vvoice = abs(Sm(:,:,4));
    Vtot=abs(X);
    
%     % Learn NMF on music only
%     Wini = rand(F,K); Hini = rand(K,T);
%     
%     [Wm_levy,Hm_levy] = levy_NMF_ME(Vmusic,Wini,Hini,Nnmf,ones(F,T),'ME');
%     Vm_levy = (Wm_levy*Hm_levy).^2;
%     
%     [Wm_is,Hm_is] = NMF(Vmusic.^2,Wini,Hini,Nnmf,0,0);
%     Vm_is = sqrt(Wm_is*Hm_is);
%     
%     [Wm_c,Hm_c] = cauchy_NMF(Vmusic,Wini,Hini,Nnmf);
%     Vm_c = Wm_c*Hm_c;
%     
%     [aux] = inexact_alm_rpca(Vmusic);
%     Vm_r = abs(aux);
%     
%     [Wm_kl,Hm_kl] = NMF(Vmusic,Wini,Hini,Nnmf,1,0);
%     Vm_kl = Wm_kl*Hm_kl;
%     
%     % Synthesis with mixture phase
%     Xm_levy = Vm_levy ./ (Vm_levy+Vvoice).* X; xm_levy = iSTFT(Xm_levy,Nfft,w,hop);
%     Xm_c = Vm_c ./ (Vm_c+Vvoice).* X; xm_c = iSTFT(Xm_c,Nfft,w,hop);
%     Xm_is = Vm_is.^2 ./ (Vm_is.^2+Vvoice.^2).* X; xm_is = iSTFT(Xm_is,Nfft,w,hop);
%     Xm_kl = Vm_kl ./ (Vm_kl+Vvoice).* X; xm_kl = iSTFT(Xm_kl,Nfft,w,hop);
%     Xm_r = Vm_r ./ (Vm_r+Vvoice).* X; xm_r = iSTFT(Xm_r,Nfft,w,hop);
%     
%     % Score
%     score(it,1) = bss_eval_sources(xtot_is',xm_is');
%     score(it,2) = bss_eval_sources(xtot_kl',xm_kl');
%     score(it,3) = bss_eval_sources(xtot_c',xm_c');
%     score(it,4) = bss_eval_sources(xtot_r',xm_r');
%     score(it,5) = bss_eval_sources(xtot_levy',xm_levy');
    
    
    % Learn NMF on music+voice
    Wini = rand(F,K); Hini = rand(K,T);
    
    [Wtot_levy,Htot_levy] = levy_NMF(Vtot,Wini,Hini,Nnmf);
    Vtot_levy = (Wtot_levy*Htot_levy).^2;   
    
    [Wtot_c,Htot_c] = cauchy_NMF(Vtot,Wini,Hini,Nnmf);
    Vtot_c = (Wtot_c*Htot_c);   
    
    [Wtot_is,Htot_is] = NMF(Vtot.^2,Wini,Hini,Nnmf,0,0);
    Vtot_is = sqrt(Wtot_is*Htot_is);
    
    [Wtot_kl,Htot_kl] = NMF(Vtot,Wini,Hini,Nnmf,1,0);
    Vtot_kl = Wtot_kl*Htot_kl;
    
    [aux] = inexact_alm_rpca(Vtot);
    Vtot_r = abs(aux);
    
    % Synthesis with mixture phase
    Xtot_levy = Vtot_levy ./ (Vtot_levy+Vvoice).* X; xtot_levy = iSTFT(Xtot_levy,Nfft,w,hop);
    Xtot_is = Vtot_is.^2 ./ (Vtot_is.^2+Vvoice.^2).* X; xtot_is = iSTFT(Xtot_is,Nfft,w,hop);
    Xtot_kl = Vtot_kl ./ (Vtot_kl+Vvoice).* X; xtot_kl = iSTFT(Xtot_kl,Nfft,w,hop);
    Xtot_c = Vtot_c ./ (Vtot_c+Vvoice).* X; xtot_c = iSTFT(Xtot_c,Nfft,w,hop);
    Xtot_r = Vtot_r ./ (Vtot_r+Vvoice).* X; xtot_r = iSTFT(Xtot_r,Nfft,w,hop);
    
    % Score
    score2(it,1) = bss_eval_sources(xtot_is',smusic);
    score2(it,2) = bss_eval_sources(xtot_kl',smusic);
    score2(it,3) = bss_eval_sources(xtot_c',smusic);
    score2(it,4) = bss_eval_sources(xtot_r',smusic);
    score2(it,5) = bss_eval_sources(xtot_levy',smusic);    
    
end

% Plot spectrograms
% figure;
% subplot(5,2,1); imagesc(ts,freq,log10(Vmusic)); axis xy; ylabel('Frequency (Hz)'); title('Music');
% subplot(5,2,2); imagesc(ts,freq,log10(Vtot)); axis xy; title('Music+Voice');
% 
% subplot(5,2,3); imagesc(ts,freq,log10(abs(Xm_is))); axis xy; title('ISNMF on Music');
% subplot(5,2,4); imagesc(ts,freq,log10(abs(Xtot_is))); axis xy; title('ISNMF on Music+Voice');
% 
% subplot(5,2,5); imagesc(ts,freq,log10(abs(Xm_kl))); axis xy; title('KLNMF on Music');
% subplot(5,2,6); imagesc(ts,freq,log10(abs(Xtot_kl))); axis xy; title('KLNMF on Music+Voice');
% 
% subplot(5,2,7); imagesc(ts,freq,log10(abs(Xm_c))); axis xy; title('Cauchy NMF on Music');
% subplot(5,2,8); imagesc(ts,freq,log10(abs(Xtot_c))); axis xy; title('Cauchy NMF on Music+Voice');
% 
% subplot(5,2,9); imagesc(ts,freq,log10(abs(Xm_levy))); axis xy; title('Lévy NMF on Music'); xlabel('Time (s)');
% subplot(5,2,10); imagesc(ts,freq,log10(abs(Xtot_levy))); axis xy; title('Lévy NMF on Music+Voice'); xlabel('Time (s)');


% Record
audiowrite(strcat(direc,'mix.ogg'),x,Fs);
audiowrite(strcat(direc,'music.ogg'),smusic,Fs);
audiowrite(strcat(direc,'music_LEVY.ogg'),xtot_levy,Fs);
audiowrite(strcat(direc,'music_IS.ogg'),xtot_is,Fs);
audiowrite(strcat(direc,'music_CAUCHY.ogg'),xtot_c,Fs);
audiowrite(strcat(direc,'music_KL.ogg'),xtot_kl,Fs);
audiowrite(strcat(direc,'music_RPCA.ogg'),xtot_r,Fs);

% 
% % Display score
figure; boxplot(score2,{'ISNMF','KLKNMF','Cauchy NMF','RPCA','Lévy NMF'});
h=ylabel('SDR (dB)'); set(h,'Fontsize',16);
ax = axis; axis([ax(1) ax(2) 5 22]);

