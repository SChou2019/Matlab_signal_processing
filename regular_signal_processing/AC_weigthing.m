clc
clear all
path = 'C:\Users\Administrator\Desktop\txt测试\test.wav';
[data,Fs] = audioread(path);
disp(Fs)
sample_points = length(data);
sam_t = sample_points/Fs;

AweightFilt = weightingFilter('A-weighting',Fs);%A计权滤波器
data_Awei = AweightFilt(data);%计权作用

CweightFilt = weightingFilter('C-weighting',Fs);
data_Cwei = CweightFilt(data);

%colormap分析,spectrogram函数
window = hann(Fs);
figure(1)
subplot(131)
%归一化显示
%spectrogram(data,window,Fs/2,Fs);
[~,F,T,P]= spectrogram(data,window,Fs/2,Fs);
%时间，频率修正

%幅值如何修正
F1 = 0:Fs/2;
n = 1:64; %64个数据快
t = (n+1)/2;%转换成时间点，以时间结束点表示
 %surf(T,F,10*log10(P),'edgecolor','none'); 
 surf(t,F1,10*log10(P),'edgecolor','none'); 
 axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
subplot(132)
%Spectrum_A = spectrogram(data_Awei,window,Fs/2,Fs);
[S,F,T,P]= spectrogram(data_Awei,window,Fs/2,Fs);
 surf(T,F,10*log10(P),'edgecolor','none'); 
 axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
subplot(133)
%Spectrum_C = spectrogram(data_Cwei,window,Fs/2,Fs);
[S,F,T,P]= spectrogram(data_Cwei,window,Fs/2,Fs);
 surf(T,F,10*log10(P),'edgecolor','none'); 
 axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');


