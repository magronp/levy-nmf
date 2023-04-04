function [bass, drums, other, vocals, accompaniment] = myfunction_rpca(mixture, fs)
warning('off','all')

 addpath(genpath('RPCA Huang'))
  [accompaniment(:,1),vocals(:,1)] = rpca_huang(mixture(:,1),fs);
  [accompaniment(:,2),vocals(:,2)] = rpca_huang(mixture(:,2),fs);
  
 mixture_length = length(mixture);
 accompaniment_length = length(accompaniment);
 vocals_length = length(vocals);
 if accompaniment_length > mixture_length
     accompaniment = accompaniment(1:mixture_length,:);
 elseif accompaniment_length < mixture_length
     accompaniment = cat(1,accompaniment,zeros(mixture_length-accompaniment_length,size(accompaniment,2)));
 end
 if vocals_length > mixture_length
     vocals = vocals(1:mixture_length,:);
 elseif accompaniment_length < mixture_length
     vocals = cat(1,vocals,zeros(mixture_length-vocals_length,size(vocals,2)));
 end
 bass = [];
 drums = [];
 other = [];

 % % Once without unvoiced lead estimation, once with it (see simm_durrieu)
% addpath(genpath('SIMM Durrieu'))
% [accompaniment,vocals,~,~,~] = simm_durrieu(mixture,fs);
% % [~,~,accompaniment,vocals,~] = simm_durrieu(mixture,fs);
% mixture_length = length(mixture);
% accompaniment_length = length(accompaniment);
% vocals_length = length(vocals);
% if accompaniment_length > mixture_length
%     accompaniment = accompaniment(1:mixture_length,:);
% elseif accompaniment_length < mixture_length
%     accompaniment = cat(1,accompaniment,zeros(mixture_length-accompaniment_length,size(accompaniment,2)));
% end
% if vocals_length > mixture_length
%     vocals = vocals(1:mixture_length,:);
% elseif accompaniment_length < mixture_length
%     vocals = cat(1,vocals,zeros(mixture_length-vocals_length,size(vocals,2)));
% end
% bass = [];
% drums = [];
% other = [];
