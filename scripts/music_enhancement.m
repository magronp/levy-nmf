%%% Enhancement of music accompaniment

clear all; clc; close all;
pkg load signal statistics;

% Random seed for reproducibility
rand("state", 1);

% Data parameters
Fs = 44100;
Ndata = 30;
dataset_path = 'data/DSD100/';
out_audio_path = 'audio_out/';

% STFT parameters
Nfft = 4096; Nw = 4096; hop = Nw/4; wtype = 'hann';

% NMF parameters
Nnmf = 200;
K = 30;

% Initialize score arrays
score = zeros(Ndata,5);
score2 = zeros(Ndata,5);

% Loop over songs
for it=1:Ndata
    clc; fprintf(' Song %d / %d \n',it,Ndata);
    
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,'Test',it,Fs,Nfft,Nw,hop,[70 80],wtype,'singing_sep');
    [F,T] = size(X);
    V = abs(Sm);
    Vtot=abs(X);
    Vmusic = V(:,:,1); Vvoice = V(:,:,2);


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
    Xtot_levy = Vtot_levy ./ (Vtot_levy+Vvoice).* X; xtot_levy = iSTFT(Xtot_levy, Nfft, hop, Nw, wtype);
    Xtot_is = Vtot_is.^2 ./ (Vtot_is.^2+Vvoice.^2).* X; xtot_is = iSTFT(Xtot_is, Nfft, hop, Nw, wtype);
    Xtot_kl = Vtot_kl ./ (Vtot_kl+Vvoice).* X; xtot_kl = iSTFT(Xtot_kl, Nfft, hop, Nw, wtype);
    Xtot_c = Vtot_c ./ (Vtot_c+Vvoice).* X; xtot_c = iSTFT(Xtot_c, Nfft, hop, Nw, wtype);
    Xtot_r = Vtot_r ./ (Vtot_r+Vvoice).* X; xtot_r = iSTFT(Xtot_r, Nfft, hop, Nw, wtype);
    
    % Score
    score2(it,1) = GetSDR(xtot_is',sm(1,:));
    score2(it,2) = GetSDR(xtot_kl',sm(1,:));
    score2(it,3) = GetSDR(xtot_c',sm(1,:));
    score2(it,4) = GetSDR(xtot_r',sm(1,:));
    score2(it,5) = GetSDR(xtot_levy',sm(1,:));    
    
    % Record audio
    rec_dir = strcat(out_audio_path,'song_',int2str(it),'/');
    mkdir(rec_dir);
    audiowrite(strcat(rec_dir,'mix.ogg'),x,Fs);
    audiowrite(strcat(rec_dir,'music.ogg'),sm(1,:),Fs);
    audiowrite(strcat(rec_dir,'music_LEVY.ogg'),xtot_levy,Fs);
    audiowrite(strcat(rec_dir,'music_IS.ogg'),xtot_is,Fs);
    audiowrite(strcat(rec_dir,'music_CAUCHY.ogg'),xtot_c,Fs);
    audiowrite(strcat(rec_dir,'music_KL.ogg'),xtot_kl,Fs);
    audiowrite(strcat(rec_dir,'music_RPCA.ogg'),xtot_r,Fs);

end

% % Display score
algos = {'ISNMF','KLKNMF','Cauchy NMF','RPCA','Lévy NMF'};
Nalgos = length(algos);

figure; boxplot(score2);
set(gca,'xtick', 1:Nalgos, 'xticklabel', algos,'fontsize',14);
h=ylabel('SDR (dB)'); set(h,'Fontsize',16);
ax = axis; axis([ax(1) ax(2) 5 22]);
