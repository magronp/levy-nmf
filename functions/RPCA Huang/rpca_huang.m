%   Huang's Robust Principal Component Analysis (RPCA) for music accompaniment/singing voice separation
%       [y1,y2] = rpca_huang(x,fs);
%
%   Input(s):
%       x: mixture [l samples, k channels]
%       fs: sampling frequency in Hz
%
%   Output(s):
%       y1: (low-rank) music accompaniment estimate [l samples, k channels]
%       y2: (sparse) singing voice estimate [l samples, k channels]
%
%   See also rpca_mask_fun, https://sites.google.com/site/singingvoiceseparationrpca/

%   Author: Zafar RAFII (zafarrafii@u.northwestern.edu)
%   Update: January 2013
%   Reference(s):
%       [1]: Po-Sen Huang and Scott Deeann Chen and Paris Smaragdis and Mark Hasegawa-Johnson. 
%            "Singing-Voice Separation from Monaural Recordings using Robust Principal Component Analysis," 
%            37th International Conference on Acoustics, Speech and Signal Processing,
%            Kyoto, Japan, March 25-30, 2012.

function [y1,y2] = rpca_huang(x,fs)

parm.lambda = 1;
parm.nFFT = 1024;
parm.windowsize = 1024;
parm.masktype = 2;                                                          %1: binary mask, 2: soft mask
parm.gain = 1;
parm.power = 1;
parm.fs = fs;

[y2,y1] = rpca_mask_fun2(x,parm);

end

%%% Modified rpca_mask_fun

function [wavoutE,wavoutA]=rpca_mask_fun2(wavinmix,parm)
    %%% parameters

    lambda = parm.lambda;
    nFFT = parm.nFFT;
    winsize = parm.windowsize;
    masktype = parm.masktype;
    gain = parm.gain;
    power = parm.power;
%     Fs= parm.fs;
%     outputname = parm.outname;

    hop = winsize/4;
    scf = 2/3;
    S = scf * stft_huang(wavinmix, nFFT ,winsize, hop);

   %%% use inexact_alm_rpca to run RPCA
    try                
        [A_mag E_mag] = inexact_alm_rpca(abs(S).^power',lambda/sqrt(max(size(S))));
        PHASE = angle(S');            
    catch err %#ok<NASGU>
        [A_mag E_mag] = inexact_alm_rpca(abs(S).^power,lambda/sqrt(max(size(S))));
        PHASE = angle(S);
    end
    
    A = A_mag.*exp(1i.*PHASE);
    E = E_mag.*exp(1i.*PHASE);

    %%% binary mask, no mask
    switch masktype                         
      case 1 % binary mask + median filter
        m= double(abs(E)> (gain*abs(A)));                  
        try  
            Emask =m.*S;
            Amask= S-Emask;
        catch err %#ok<NASGU>
            Emask =m.*S';
            Amask= S'-Emask;
        end        
      case 2 % no mask
        Emask=E;
        Amask=A;
      otherwise 
          fprintf('masktype error\n');
    end

    %%% do istft
    try 
        wavoutE = istft_huang(Emask', nFFT ,winsize, hop)';   
        wavoutA = istft_huang(Amask', nFFT ,winsize, hop)';
    catch err %#ok<NASGU>
        wavoutE = istft_huang(Emask, nFFT ,winsize, hop)';   
        wavoutA = istft_huang(Amask, nFFT ,winsize, hop)';
    end

    wavoutE=wavoutE/max(abs(wavoutE));
%     wavwrite(wavoutE,Fs,[outputname,'_E']);

    wavoutA=wavoutA/max(abs(wavoutA));
%     wavwrite(wavoutA,Fs,[outputname,'_A']);

%     %% evaluate
%     if length(wavoutA)==length(wavinA)
% 
%         sep = [wavoutA , wavoutE]';
%         orig = [wavinA , wavinE]';
% 
%         for i = 1:size( sep, 1)
%                [e1,e2,e3] = bss_decomp_gain( sep(i,:), i, orig);
%                [sdr(i),sir(i),sar(i)] = bss_crit( e1, e2, e3);
%         end
%     else
%         minlength=min( length(wavoutE), length(wavinE) );
% 
%         sep = [wavoutA(1:minlength) , wavoutE(1:minlength)]';
%         orig = [wavinA(1:minlength) , wavinE(1:minlength)]';
% 
%         for i = 1:size( sep, 1)
%                [e1,e2,e3] = bss_decomp_gain( sep(i,:), i, orig);
%                [sdr(i),sir(i),sar(i)] = bss_crit( e1, e2, e3);
%         end
%     end
% 
%     Parms.SDR=sdr(2);
%     Parms.SIR=sir(2);
%     Parms.SAR=sar(2);

end
