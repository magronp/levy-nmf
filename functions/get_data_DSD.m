function [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,data_split,ind,Fs,Nfft,Nw,hop,t_chunk,wtype,task)

Fs_old=44100; % original sample rate

if nargin<10
    task = 'all_sources';
end

if nargin<9
    wtype = 'hann';
end

% Specify a time chunk in seconds
if nargin<8
    t_chunk = [70 80];
end

% Identify and skip empty songs
Nsongs = 50;
datavec = 1:Nsongs;
switch data_split
    case 'Dev'
        dataNaN=[1 18 29 38 49];
    case 'Test'
        dataNaN=[6 34 36 40];
end
datavec(dataNaN)=[];
datavec = datavec(1:min(Nsongs,length(datavec)));
num_piece = datavec(ind);

%%% Get files path
dir_path = strcat(dataset_path,'Sources/',data_split);
aux = dir(dir_path);
song_name = aux(num_piece+2).name;
list_instr = {'bass','drums','other','vocals'};
J = length(list_instr);
L = strcat(dir_path,'/',song_name,'/',list_instr,'.wav');

%%% Read time-domain signals

t_beg = t_chunk(1); t_end = t_chunk(2);
sig_length = ceil(((t_end-t_beg)+1/Fs_old)*Fs);
s_aux=zeros(J,sig_length);
for j = 1:J
    aux = audioread(L{j},Fs_old*[t_beg t_end]);   % read source
    aux=mean(aux,2)';                             % average channels
    s_aux(j,:) = resample(aux,Fs,Fs_old);         % resample
end

%% If singing voice separation, form new sources
if strcmp(task,'singing_sep')
    J=2;
    s_aux = [s_aux(1,:)+s_aux(2,:)+s_aux(3,:) ; s_aux(4,:)];
end


%%% STFT and iSTFT to ensure proper signal length
Sm = STFT(s_aux,Nfft,hop,Nw,wtype);
sm = iSTFT(Sm,Nfft,hop,Nw,wtype);
x = sum(sm,1); X = sum(Sm,3);

ts = (0:size(X,2)-1)*hop / Fs;
freq = (1:Nfft/2)*Fs/Nfft;


end