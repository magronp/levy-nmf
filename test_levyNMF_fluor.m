clear all; close all; clc;
nNMF = 50;
Nal=3;

data_dir = 'data/';

% Load fluorescence data
load(strcat(data_dir,'fluor/ble.mat'));

X=Grain_ble_EUSIPCO_2004.Data; X=X';
[F,T] = size(X);
K=3; N=20;
labs={'Acide félurique lié','Acide félurique libre','Acide paracoumarique'}; longond = 'Longueur d''onde (nm)';
labs={'Bound feluric acid','Free feluric acid','p-coumaric acid'}; longond = 'Wavelength (nm)';

freq = Grain_ble_EUSIPCO_2004.WvL;
pos_xy=Grain_ble_EUSIPCO_2004.positions_xy;
pos_x=Grain_ble_EUSIPCO_2004.positions_x;
pos_y=Grain_ble_EUSIPCO_2004.positions_y;

% Oracle: load pure spectra and learn activations
or = load(strcat(data_dir,'fluor/spectre_pur.mat'));
Wo = or.data;

Wini=rand(F,K); Hini=rand(K,T);
[~,Ho] = NMF(X.^2,Wo,Hini,nNMF,0,0,0,ones(F,T),0);
Xo = wiener(X,Wo,Ho);

% NMF (euclidean, KL and Lévy)
[We,He] = NMF(X,Wo,Hini,nNMF,2,0,0);
[Wk,Hk] = NMF(X,Wo,Hini,nNMF,1,0,0);
[Wl,Hl] = levy_NMF(X,Wo,Hini,nNMF);

% Sort estimated with oracle
[We,He] = sort_distcorr(Wo,We,He);
[Wk,Hk] = sort_distcorr(Wo,Wk,Hk);
[Wl,Hl] = sort_distcorr(Wo,Wl,Hl);


% Source separation (generalized Wiener)
Xe=zeros(F,T,K,Nal);
Xe(:,:,:,1) = wiener(X,We,He); Xe(:,:,:,2) = wiener(X,Wk,Hk); Xe(:,:,:,3) = wiener(X,Wl,Hl);


% Estimation error
err_kl = zeros(K,Nal); cor=zeros(K,Nal);
for al=1:Nal
    for k=1:K
        err_kl(k,al) = beta_div(squeeze(Xo(:,:,k)),squeeze(Xe(:,:,k,al)),1); %KL div
        cor(k,al) = distcorr(squeeze(Xo(:,:,k)),squeeze(Xe(:,:,k,al)));
    end
end

% Plot KL
figure;
bar(log(err_kl));
set(gca,'XTickLabel',labs);
ha=ylabel('$\log (d_{KL})$');set(ha,'FontSize',16,'interpreter','latex');
ha=legend('Euc','KL','Lévy'); set(ha,'FontSize',14);

% Plot corr
figure;
bar(cor);
set(gca,'XTickLabel',labs);
ha=ylabel('Corrélation');set(ha,'FontSize',16,'interpreter','latex');
ha=legend('EuNMF','KLNMF','Lévy NMF'); set(ha,'FontSize',14);



% Concentrations of components (Oracle case)
concentr = zeros(N,N,K);
figure; colormap(gray);
for k=1:K
    concentr(:,:,k) = reshape(Hl(k,:),[N N]);
    %subplot(1,3,k); imagesc(pos_x,pos_y,concentr(:,:,k)); axis xy; xlabel('x (nm)'); ylabel('y (nm)'); title(labs{k});
    subplot(1,3,k); imagesc(concentr(:,:,k)); axis xy; title(labs{k}); set(gca,'Xtick',[],'Ytick',[]);
end

% Plot bases
We = We ./ repmat(max(We)+eps,[F 1]);
Wl = Wl ./ repmat(max(Wl)+eps,[F 1]);
Wk = Wk ./ repmat(max(Wk)+eps,[F 1]);
Wo = Wo ./ repmat(max(Wo)+eps,[F 1]);

figure;
for k=1:3
    subplot(1,3,k);
    plot(freq,Wo(:,k),'k'); hold on;
    plot(freq,We(:,k),'b--'); hold on;
    plot(freq,Wk(:,k),'r-.'); hold on;
    plot(freq,Wl(:,k),'m.-'); 
    axis([freq(1) freq(end) 0 1]); ha=xlabel(longond); set(ha,'FontSize',16); title(labs{k});
end
ha=legend('Original','Euc','KL','Lévy'); set(ha,'FontSize',14);

