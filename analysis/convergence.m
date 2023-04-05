%%% Convergence of Lévy NMF algorithms (heuristic, MM, and ME)

clear all; clc; close all;
pkg load signal statistics;

% Random seed for reproducibility
rand("state", 1);

% Data parameters
Fs = 8000;
dataset_path = 'data/guitar/';

% STFT parameters
Nfft = 4096; Nw = 4096; hop = Nw/4; wtype = 'hann';

% NMF
Nnmf = 100;
K = 5;

% Load the data
[x,X,F,T,ts,freq] = get_data_guitar(dataset_path,1,Fs,Nfft,hop,Nw,wtype);
V = abs(X);

% Init NMF
Wini = rand(F,K); Hini = rand(K,T);
mask = ones(F,T);

% Levy NMF: Heuristic, MM and ME algorithms
[~,~,err_h] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'H');
[~,~,err_mm] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'MM');
[~,~,err_me1] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','line');
[~,~,err_me2] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','DL');
[~,~,err_me3] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','newton');
[~,~,err_me4] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','dicho');

% Compare heuristic, MM and ME
semilogy(0:Nnmf,err_h,'k.-',0:Nnmf,err_mm,'b*-',0:Nnmf,err_me1,'r-');
ha=legend('Heuristic','MM','ME'); set(ha,'FontSize',14);
ha=xlabel('Iterations'); set(ha,'FontSize',16);
ha=ylabel('Error'); set(ha,'FontSize',16);

% Compare ME methods
figure;
semilogy(1:Nnmf,err_me1(2:end),'k.-',1:Nnmf,err_me2(2:end),'b*-',1:Nnmf,err_me3(2:end),'ro-',1:Nnmf,err_me4(2:end),'md-');
ha=legend('line','DL','newton','dicho'); set(ha,'FontSize',14);
ha=xlabel('Iterations'); set(ha,'FontSize',16);
ha=ylabel('Error'); set(ha,'FontSize',16);
