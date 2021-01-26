
pink_wave = create_pink_noise(51200,5,-10);
sound(pink_wave,51200)

function out = create_pink_noise(Fs, Sec, Amp)
% Creates a pink noise signal and saves it as a wav file
%
% Usage: create_noise(Fs, Sec, Amp);
%
%        Fs is the desired sampling rate
%        Sec is the duration of the signal in seconds
%        Amp is the amplitude in dB of the signal (0dB to -144dB)
%
% Author: sparafucile17 06/14/02
%error trapping
if((Amp > 0) || (Amp < -144))
error('Amplitude is not within the range of 0dB to -144dB');
end
%Create Whitenoise
white_noise = randn((Fs*Sec)+1,1);
%Apply weighted sum of first order filters to approximate a -10dB/decade
%filter.  This is Paul Kellet's "refined" method (a.k.a instrumentation
%grade)  It is accurate to within +/-0.05dB above 9.2Hz
b=zeros(7,1);
for i=1:((Fs*Sec)+1)
b(1) = 0.99886 * b(1) + white_noise(i) * 0.0555179;
b(2) = 0.99332 * b(2) + white_noise(i) * 0.0750759;
b(3) = 0.96900 * b(3) + white_noise(i) * 0.1538520;
b(4) = 0.86650 * b(4) + white_noise(i) * 0.3104856;
b(5) = 0.55000 * b(5) + white_noise(i) * 0.5329522;
b(6) = -0.7616 * b(6) - white_noise(i) * 0.0168980;
pink_noise(i) = b(1) + b(2) + b(3) + b(4) + b(5) + b(6) + b(7) + white_noise(i) * 0.5362;
b(7) = white_noise(i) * 0.115926; 
end
%Normalize to +/- 1
if(abs(min(pink_noise)) > max(pink_noise))
    pink_noise = pink_noise / abs(min(pink_noise));
else
    pink_noise = pink_noise / max(pink_noise);
end
%Normalize to prevent positive saturation (We can't represent +1.0)
pink_noise = pink_noise /abs(((2^31)-1)/(2^31));
%Scale signal to match desired level
pink_noise = pink_noise * 10^(Amp/20);
%Output noise signal
out = pink_noise(1:end-1);
end
