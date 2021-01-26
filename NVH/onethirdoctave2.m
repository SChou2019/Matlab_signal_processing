clc;clear all
%[x, Fs] = audioread('C:\Users\Administrator\Desktop\export\Noise_PinkNoise_2_115_53_51.wav'); 
path = "C:\Users\Administrator\Desktop\export\Noise_PinkNoise_2_115_52_57.mat";
x = load(path);
x = x.Noise_PinkNoise_2_115_52_57;
Fs = 32768;
Fc =  [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 ...
500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 ... 
6300 8000 10000 12500 16000 20000];
Fc = Fc(Fc<Fs/2);
center_num = length(Fc);
fl = Fc*2^(-1/6);%lower frequency
fu = Fc*2^(1/6);%higher frequency
fu(fu>Fs/2) = Fs/2-1;
data = x(1:Fs);%取1s数据
len = length(data);
fft_amp = abs(fft(data))*2/len;
fft_amp = fft_amp(1:floor(len/2));
fft_amp_dB = 20*log10(fft_amp/2e-5);
delta_f = Fs/len;
f = 1:delta_f:Fs/2;
oct_amp = zeros(center_num,1);
for i = 1:center_num
    lower = fl(i);
    upper = fu(i);
%     for j = 1:length(f)
%         if f(j)>lower && f(j)<=upper
%             oct_amp(i) = oct_amp(i) + fft_amp(j)^2;
%         end
%     end
    index = find(f>lower & f<upper);
    %index = find(f(i)>lower & f(i)<upper);%继续排查该函数存在的问题
    for j = 1:length(index)
        oct_amp(i) = oct_amp(i) + fft_amp(index(j))^2;
    end
end
oct_dB = 10*log10(oct_amp/2e-5);

figure(1)
plot(oct_dB)
% xlim([0 30])