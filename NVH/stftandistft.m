clear, clc, close all

% load an audio signal
[x, fs] = audioread('track.wav'); 
x = x(:, 1);                                  

% signal parameters
xlen = length(x);                   
t = (0:xlen-1)/fs;                  

% define the analysis and synthesis parameters
wlen = 1024;
hop = wlen/8;
nfft = 4*wlen;

% generate analysis and synthesis windows
anal_win = blackmanharris(wlen, 'periodic');
synth_win = hamming(wlen, 'periodic');

% perform time-frequency analysis and resynthesis of the signal
[STFT, ~, ~] = stft(x, anal_win, hop, nfft, fs);
[x_istft, t_istft] = istft(STFT, anal_win, synth_win, hop, nfft, fs);

% plot the original signal
figure(1)
plot(t, x, 'b')
grid on
xlim([0 max(t)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Signal amplitude')
title('Original and reconstructed signal')

% plot the resynthesized signal 
hold on
plot(t_istft, x_istft, '-.r')
legend('Original signal', 'Reconstructed signal')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Inverse Short-Time Fourier Transform        %
%               with MATLAB Implementation             %
%                                                      %
% Author: Ph.D. Eng. Hristo Zhivomirov        12/26/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, t] = istft(STFT, awin, swin, hop, nfft, fs)

% function: [x, t] = istft(STFT, awin, swin, hop, nfft, fs)
%
% Input:
% stft - STFT-matrix (only unique points, time
%        across columns, frequency across rows)
% awin - analysis window function
% swin - synthesis window function
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
%
% Output:
% x - signal in the time domain
% t - time vector, s

% signal length estimation and preallocation
L = size(STFT, 2);          % determine the number of signal frames
wlen = length(swin);        % determine the length of the synthesis window
xlen = wlen + (L-1)*hop;    % estimate the length of the signal vector
x = zeros(1, xlen);         % preallocate the signal vector

% reconstruction of the whole spectrum
if rem(nfft, 2)             
    % odd nfft excludes Nyquist point
    X = [STFT; conj(flipud(STFT(2:end, :)))];
else                        
    % even nfft includes Nyquist point
    X = [STFT; conj(flipud(STFT(2:end-1, :)))];
end

% columnwise IFFT on the STFT-matrix
xw = real(ifft(X));
xw = xw(1:wlen, :);

% Weighted-OLA
for l = 1:L
    x(1+(l-1)*hop : wlen+(l-1)*hop) = x(1+(l-1)*hop : wlen+(l-1)*hop) + ...
                                      (xw(:, l).*swin)';
end

% scaling of the signal
W0 = sum(awin.*swin);                  
x = x.*hop/W0;                      

% generation of the time vector
t = (0:xlen-1)/fs;                 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Short-Time Fourier Transform            %
%               with MATLAB Implementation             %
%                                                      %
% Author: Ph.D. Eng. Hristo Zhivomirov        12/21/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [STFT, f, t] = stft(x, win, hop, nfft, fs)

% function: [STFT, f, t] = stft(x, win, hop, nfft, fs)
%
% Input:
% x - signal in the time domain
% win - analysis window function
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
%
% Output:
% STFT - STFT-matrix (only unique points, time 
%        across columns, frequency across rows)
% f - frequency vector, Hz
% t - time vector, s

% representation of the signal as column-vector
x = x(:);

% determination of the signal length 
xlen = length(x);

% determination of the window length
wlen = length(win);

% stft matrix size estimation and preallocation
NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
STFT = zeros(NUP, L);       % preallocate the stft matrix

% STFT (via time-localized FFT)
for l = 0:L-1
    % windowing
    xw = x(1+l*hop : wlen+l*hop).*win;
    
    % FFT
    X = fft(xw, nfft);
    
    % update of the stft matrix
    STFT(:, 1+l) = X(1:NUP);
end

% calculation of the time and frequency vectors
t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
f = (0:NUP-1)*fs/nfft;

end