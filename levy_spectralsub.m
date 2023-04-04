clear all;

% Get data
[s,Fs_old] = audioread('datasets/guitar/pieces/lick1.wav'); s = s';
Fs = 8000;
s = resample(s,Fs,Fs_old);
s = s(1.5*Fs:round(9.3*Fs));
N = length(s);

% Add noise
SNR = 0;
Ex = real(s*s') / N;
b = randn(1,N);
Eb = real(b*b') / N;
b = b * sqrt(Ex/Eb) * 10^(-SNR/20);
x = s + b;

% STFT
Nfft = 1000; hop = Nfft/4;
X = STFT(x,Nfft,hop);  S = STFT(s,Nfft,hop); B = STFT(b,Nfft,hop);
[F,T] = size(X);

%Spectral subtraction
alpha = 2; beta = 7;
noise_time_av = mean(abs(B),2);
est_noise = repmat(noise_time_av,[1 T]);

Vhat = max(0,abs(X).^alpha - beta * est_noise.^alpha).^(1/alpha);

% Levy NMF for musical noise removal
nNMF = 200; K = 30;
Wini = rand(F,K); Hini = rand(K,T);

[Wk,Hk] = NMF(Vhat,Wini,Hini,nNMF,1,0,0,ones(F,T));
[Wl,Hl] = levy_NMF(Vhat,Wk,Hini,nNMF);

Vlevy = (Wl*Hl);
Vkull = (Wk*Hk).^2;

% Plot spectro
subplot(2,3,1); imagesc((abs(S))); axis xy;
subplot(2,3,2); imagesc((abs(X))); axis xy;
subplot(2,3,3); imagesc((Vhat)); axis xy;
subplot(2,3,4); imagesc((Vkull)); axis xy;
subplot(2,3,5); imagesc((Vlevy)); axis xy;

% Synthesis
Shat = Vhat .* exp(1i * angle(X)); shat = iSTFT(Shat,Nfft,hop);
Slevy = Vlevy .* exp(1i * angle(X)); slevy = iSTFT(Slevy,Nfft,hop);
Skull = Vkull .* exp(1i * angle(X)); skull = iSTFT(Skull,Nfft,hop);


%soundsc(shat,Fs);
soundsc(slevy,Fs);






