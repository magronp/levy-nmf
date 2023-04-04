clear all; clc; close all;
Fs = 44100; Nfft = 4096; Nw = 4096; hop = Nw/4;
source_type='PIANO_FIFTH';
Nnmf=200;
Ninit=30;
Ndata=30;
Nal=4;

errorsc_std = zeros(Ndata,4);

for it=1:Ndata

    midi_pitch=[30 37]+it;
    % Generate sources
    num_piece = it; gen_sources_time;
    gen_sources_TF;
    %gen_sources_TF_inst; K=40;
    V = abs(X);
    err = zeros(Ninit,4);
    err_kl = zeros(Ninit,4);
    errsc = zeros(Ninit,4);
    for init=1:Ninit
        clc; fprintf('data %d / %d \n init %d / %d \n',it,Ndata,init,Ninit);

        % New initialization
        Wini = rand(F,K); Hini = rand(K,T);

        % NMFs
        [Wis,His] = NMF(V.^2,Wini,Hini,Nnmf,0,0);
        [Wkl,Hkl] = NMF(V,Wini,Hini,Nnmf,1,0);
        [Wc,Hc] = cauchy_NMF(V,Wini,Hini,Nnmf);
        [Wl,Hl] = levy_NMF_ME(V,Wini,Hini,Nnmf,ones(F,T),'ME');
        
        % Wiener
        Xis = wiener(V,Wis,His);
        Xkl = wiener(V,Wkl,Hkl);
        Xc = wiener(V,Wc,Hc);
        Xl = wiener(V,Wl,Hl);
        Xe = zeros(F,T,K,Nal); Xe(:,:,:,1)=Xis;Xe(:,:,:,2)=Xkl;Xe(:,:,:,3)=Xc;Xe(:,:,:,4)=Xl;
        
        % Error
        aux = zeros(K,4);
        for k=1:K
            for al=1:Nal
                aux(k,al) = norm(abs(Sm(:,:,k))-squeeze(Xe(:,:,k,al)));
            end
        end
        errsc(init,:) = mean(aux,1);
    end
    
    errorsc_std(it,:) = std(errsc);

end

% Display score
labs = {'IS','KL','Cauchy','Lévy'};
bar(mean(errorsc_std));
set(gca,'XTickLabel',labs);
ylabel('Écart-type moyen');



