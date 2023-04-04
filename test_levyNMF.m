clear all; clc; close all;
Fs = 44100; Nfft = 512*8; Nw = 512*8; hop = Nw/4;
Nnmf=100;

source_type='SYNTH_NO_OL';
gen_sources_time; gen_sources_TF;
%gen_sources_TF_inst;

% Nonnegative data
V = abs(X);
Wini = rand(F,K); Hini = rand(K,T);
mask = ones(F,T);

% Levy NMF: Heuristic, MM and ME algorithms
[~,~,err_h] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'H');
[~,~,err_mm] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'MM');

[~,~,err_me1] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','line');
[~,~,err_me2] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','DL');
[~,~,err_me3] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','newton');
[~,~,err_me4] = levy_NMF_ME(V,Wini,Hini,Nnmf,mask,'ME','dicho');


% Compare H, MM and ME
semilogy(0:Nnmf,err_h,'k.-',0:Nnmf,err_mm,'b*-',0:Nnmf,err_me1,'r-');
ha=legend('Naïve','MM','ME'); set(ha,'FontSize',14);
ha=xlabel('Itérations'); set(ha,'FontSize',16);
ha=ylabel('Erreur'); set(ha,'FontSize',16);


% Compare ME methods
figure;
semilogy(1:Nnmf,err_me1(2:end),'k.-',1:Nnmf,err_me2(2:end),'b*-',1:Nnmf,err_me3(2:end),'ro-',1:Nnmf,err_me4(2:end),'md-');
ha=legend('line','DL','newton','dicho'); set(ha,'FontSize',14);
ha=xlabel('Itérations'); set(ha,'FontSize',16);
ha=ylabel('Erreur'); set(ha,'FontSize',16);


