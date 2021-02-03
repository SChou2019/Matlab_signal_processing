%ʵ��Ҫ��ADPCM�����ʵ��
clear all;
clc
close all;
[x,fs,numbits]= wavread('C6_3_y.wav'); 

sign_bit=2;                                     %��λADPCM�㷨
ss=adpcm_encoder(x,sign_bit);
yy=adpcm_decoder(ss,sign_bit)';

nq=sum((x-yy).*(x-yy))/length(x);
sq=mean(yy.^2);
snr=(sq/nq);
t=(1:length(x))/fs;
subplot(211)
plot(t,x/max(abs(x)))
axis tight
title('(a)����ǰ����')
xlabel('ʱ��/s')
ylabel('����')
subplot(212)
plot(t,yy/max(abs(yy)))
axis tight
title('(b)���������')
xlabel('ʱ��/s')
ylabel('����')
snrq=10*log10(mean(snr))


