%%% Test the capacity of Lévy NMF to model impulsive noise

clear all; clc; close all;
pkg load signal statistics;

% Random seed for reproducibility
rand("state", 1);

% Data parameters
F=50; T=50; K=5;
Ndata = 200;

% NMF parameters
Nnmf=200;
alpha = 0.1:0.01:0.5;
Nalpha=length(alpha);

% Initialize loss arrays
alpha_disp = zeros(Nalpha,5,Ndata);
kl_div = zeros(Nalpha,5,Ndata);

% Loop over data samples
for it=1:Ndata
    
    % Generate sources
    W=abs(randn(F,K)).^4; H=abs(randn(K,T)).^4;
    WH = W*H;
        
    for a=1:Nalpha    
        clc; fprintf(' data %d / %d \n alpha %d / %d \n',it,Ndata,a,Nalpha);

        % Generate alpha-stable data
        al =alpha(a);
        sigma = (WH).^(1/al);
        V = abs(stablernd(al,1,sigma,0,F,T));

        % Estimate the scale parameters
        Wini = rand(F,K); Hini = rand(K,T);
        [Wl,Hl] = levy_NMF(V,Wini,Hini,Nnmf);
        sigma_levy = (Wl*Hl).^2;

        [Wis,His] = NMF(V.^2,Wini,Hini,Nnmf,0,0);
        sigma_is = (Wis*His).^0.5;

        [Wkl,Hkl] = NMF(V,Wini,Hini,Nnmf,1,0);
        sigma_kl = Wkl*Hkl;

        [Wc,Hc,errc] = cauchy_NMF(V,Wini,Hini,Nnmf);
        sigma_c = Wc*Hc;

        [aux] = inexact_alm_rpca(V,1/sqrt(F),0,Nnmf);
        sigma_r=abs(aux);

        % Error: alpha-dispersion and KL divergence
        alpha_disp(a,1,it) = sum(abs(sigma(:)-sigma_is(:)).^(1/al));  kl_div(a,1,it) = beta_div(sigma,sigma_is,1);
        alpha_disp(a,2,it) = sum(abs(sigma(:)-sigma_kl(:)).^(1/al));  kl_div(a,2,it) = beta_div(sigma,sigma_kl,1);
        alpha_disp(a,3,it) = sum(abs(sigma(:)-sigma_c(:)).^(1/al));  kl_div(a,3,it) = beta_div(sigma,sigma_c,1);
        alpha_disp(a,4,it) = sum(abs(sigma(:)-sigma_r(:)).^(1/al));  kl_div(a,4,it) = beta_div(sigma,sigma_r,1);
        alpha_disp(a,5,it) = sum(abs(sigma(:)-sigma_levy(:)).^(1/al));  kl_div(a,5,it) = beta_div(sigma,sigma_levy,1);

    end
end

alpha_disp_av = mean(alpha_disp,3);
kl_div_av = mean(kl_div,3);

figure;
subplot(2,1,1);
plot(alpha,log(alpha_disp_av(:,1)),'g.-',alpha,log(alpha_disp_av(:,2)),'ro-',alpha,log(alpha_disp_av(:,3)),'kd-',alpha,log(alpha_disp_av(:,4)),'m+-',alpha,log(alpha_disp_av(:,5)),'b*-');
h=xlabel('\alpha'); set(h,'fontsize',16); h=ylabel('$\log (L_\alpha)$'); set(h,'fontsize',16,'interpreter','latex');
subplot(2,1,2);
plot(alpha,log(kl_div_av(:,1)),'g.-',alpha,log(kl_div_av(:,2)),'ro-',alpha,log(kl_div_av(:,3)),'kd-',alpha,log(kl_div_av(:,4)),'m+-',alpha,log(kl_div_av(:,5)),'b*-');
h=xlabel('\alpha'); set(h,'fontsize',16); h=ylabel('$\log (KL)$'); set(h,'fontsize',16,'interpreter','latex');
h=legend('ISNMF','KLNMF','Cauchy NMF','RPCA','Lévy NMF'); set(h,'fontsize',14);
