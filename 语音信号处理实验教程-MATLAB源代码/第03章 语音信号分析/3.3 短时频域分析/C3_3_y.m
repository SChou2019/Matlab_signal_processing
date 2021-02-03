%ʵ��Ҫ��������ڶ�ʱ����Ҷ����������ͼ��ʾ
clear all;
clc; 
close all;

[x,fs]=wavread('C3_3_y.wav');       % ���������ļ�
wlen=256;
nfft=wlen;
win=hanning(wlen);
inc=128;          % ����֡����֡��

y=STFFT(x,win,nfft,inc);        %���ʱ����Ҷ�任

fn=size(y,2);                           %֡��

freq=(0:wlen/2)*fs/wlen;                % ����FFT���Ƶ�ʿ̶�

frameTime=FrameTimeC(fn,wlen,inc,fs); % ����ÿ֡��Ӧ��ʱ��
imagesc(frameTime,freq,20*log10(abs(y)+eps)); % ����Y��ͼ��  
axis xy; ylabel('Ƶ��/Hz');xlabel('ʱ��/s');
title('������ͼ');
colormap(jet)
