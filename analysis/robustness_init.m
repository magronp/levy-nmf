%%% Robustness to initialization

clear all; clc; close all;
pkg load signal statistics;

% Random seed for reproducibility
rand("state", 1);

% Data parameters
Fs = 8000;
Ndata = 6;
dataset_path = 'data/guitar/';

% STFT parameters
Nfft = 4096; Nw = 4096; hop = Nw/4; wtype = 'hann';

% NMF parameters
Nnmf = 200;
Ninit = 30;
Nal = 4;
K = 30;

% Initialize loss array
reco_loss = zeros(Ndata,Ninit,4);

% Loop over data samples
for m=1:Ndata

    % Load the data
    [x,X,F,T,ts,freq] = get_data_guitar(dataset_path,m,Fs,Nfft,hop,Nw,wtype);
    V = abs(X);
    
    % Initialize arrays
    errsc = zeros(Ninit,4);
    
    % Loop over initializations
    for init=1:Ninit
        clc; fprintf('data %d / %d \n init %d / %d \n',m,Ndata,init,Ninit);

        % New initialization
        Wini = rand(F,K); Hini = rand(K,T);

        % NMFs
        [Wis,His] = NMF(V.^2,Wini,Hini,Nnmf,0,0);
        [Wkl,Hkl] = NMF(V,Wini,Hini,Nnmf,1,0);
        [Wc,Hc] = cauchy_NMF(V,Wini,Hini,Nnmf);
        [Wl,Hl] = levy_NMF_ME(V,Wini,Hini,Nnmf,ones(F,T),'ME');
        
        % Reconstruction error
        reco_loss(m,init,1) = norm(V - Wis*His);
        reco_loss(m,init,2) = norm(V - Wkl*Hkl);
        reco_loss(m,init,3) = norm(V - Wc*Hc);
        reco_loss(m,init,4) = norm(V - Wl*Hl);
    end
end

% Get the std over initializations
reco_loss_std = squeeze(std(reco_loss, [], 2));

% Display average std
labs = {'IS','KL','Cauchy','Lévy'};
bar(mean(reco_loss_std));
set(gca,'XTickLabel',labs);
ylabel('Average STD');
