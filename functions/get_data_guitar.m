function [x,X,F,T,ts,freq] = get_data_guitar(dataset_path,num_piece,Fs,Nfft,hop,Nw,wtype)

%%% Read audio data
L = strcat(dataset_path,'lick',int2str(num_piece),'.wav');
[s,Fs_old] = audioread(L);
s=mean(s,2);
xaux = resample(s,Fs,Fs_old)';

%%% STFT
X = STFT(xaux, Nfft, hop, Nw, wtype);
x = iSTFT(X, Nfft, hop, Nw, wtype);

[F,T] = size(X);
ts = (0:T-1)*hop / Fs;
freq = (1:Nfft/2)*Fs/Nfft;

end
