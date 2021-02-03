%ʵ��Ҫ��һ������ģ��
clc
clear all
close all
T1=5;                       %5ms
T2=3;                       %3ms
fs=8;                       %8kHz
N1=T1*fs;
N2=T2*fs;
n1=1:N1-1;
n2=N1:N1+N2;
g=zeros(1,20*fs);
t=(1:20*fs)/fs;
g(1:N1-1)=0.5*(1-cos(pi*n1/N1));
g(N1:N1+N2)=cos(pi*(n2-N1)/(2*N2));
plot(t,g)
ylim([-0.4,1.2])
line([0 20],[0 0])
line([5 5],[0 1],'LineStyle','--')
xlabel('ʱ��/ms')
ylabel('g(n)')
title('(a)ʱ����')

figure
N=length(g);
f=(0:N-1)*fs/N;
fftg=fft(g);
disg=20*log10(abs(fftg));
plot(f(1:N/2+1),disg(1:N/2+1))
line([0 fs/2],[0 0])
xlabel('Ƶ��/kHz')
ylabel('20lgG')
title('(b)Ƶ�ײ���')