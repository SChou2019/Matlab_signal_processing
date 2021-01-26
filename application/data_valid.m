%������֤
clc
clear all
fs = 500;
w = 2*pi*200;
t = 0:1/fs:20;
%data_example = cos(fs*t+pi/6);
data_example = 2*cos(w*t+pi/6);
data = data_example(1:1024);

AweightFilt = weightingFilter('A-weighting',fs);%A��Ȩ�˲���
data_Awei = AweightFilt(data);%��Ȩ����

original = fft(data);
original = abs(original)*2/length(data);

Awei = fft(data_Awei);
Awei = abs(Awei)*2/length(data_Awei);

f = (1:length(data)/2)*fs/length(data);
figure(1)
subplot(121)
plot(f,original(1:length(data)/2))
subplot(122)
plot(f,Awei(1:length(data)/2))

% figure
% plot(data_example)
window = hann(fs);
figure(1)
subplot(131)
%��һ����ʾ
spectrogram(data_example,window,fs/2,fs);
[~,F,T,P]= spectrogram(data_example,window,fs/2,fs);
%ʱ�䣬Ƶ������
F1 = 0:Fs/2;
n = 1:64; %64�����ݿ�
t = (n+1)/2;%ת����ʱ��㣬��ʱ��������ʾ
 %surf(T,F,10*log10(P),'edgecolor','none'); 
 surf(t,F1,10*log10(P),'edgecolor','none'); 
 axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
subplot(132)