% ʵ��Ҫ������Ԥ���ع��ܲ���
clc
clear all
close all
[s,fs]=wavread('C2_5_y_3.wav');
e=s(2000:2225);                                 %��ȡһ�ν��з��������׿����仯
un=filter([1,-0.95],1,e);                       %Ԥ�����ź�b=[1,-0.95];

%ԭʼ�ź�Ƶ��
N=512;
pinlv=(0:1:N/2-1)*fs/N;
x=fft(e,N);
r1=abs(x);
t1=20*log10(r1);
signal=t1(1:N/2);

%Ԥ�����ź�Ƶ��
[h1,w1]=freqz([1,-0.95],1,256,fs);
pha=angle(h1);
H1=abs(h1);
r2=r1(1:N/2);
u=r2.*h1;
u2=abs(u);
signalPre=20*log10(u2);

figure(1);
subplot(211)
plot(e,'b*-')
ylim([-0.4,1])
hold on
plot(real(un),'ro-')
legend('ԭʼ�����ź�','Ԥ���غ�������ź�')
title('ԭʼ�����źź�Ԥ���غ�������ź�');
xlabel('������');ylabel('����');
subplot(212);
plot(pinlv,signal,'g+-')
hold on
plot(pinlv,signalPre,'kx-')
legend('ԭʼ�����ź�Ƶ��','Ԥ���غ�������ź�Ƶ��')
title('Ԥ����ǰ��������ź�Ƶ��');
xlabel('Ƶ��');ylabel('����/dB');

