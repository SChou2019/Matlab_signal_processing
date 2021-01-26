clc
clear all
path = 'C:\Users\Administrator\Desktop\txt����\test.wav';
[data,Fs] = audioread(path);
disp(Fs)
sample_points = length(data);
sam_t = sample_points/Fs;

AweightFilt = weightingFilter('A-weighting',Fs);%A��Ȩ�˲���
data_Awei = AweightFilt(data);%��Ȩ����

CweightFilt = weightingFilter('C-weighting',Fs);
data_Cwei = CweightFilt(data);

%colormap����,spectrogram����
window = hann(Fs);
figure(1)
subplot(131)
%��һ����ʾ
%spectrogram(data,window,Fs/2,Fs);
[~,F,T,P]= spectrogram(data,window,Fs/2,Fs);
%ʱ�䣬Ƶ������

%��ֵ�������
F1 = 0:Fs/2;
n = 1:64; %64�����ݿ�
t = (n+1)/2;%ת����ʱ��㣬��ʱ��������ʾ
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


