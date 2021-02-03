%ʵ��Ҫ��һ�����׼�������ʾ
clear all; clc; close all;
[y,fs]=wavread('C3_4_y_1.wav');
y=y(1:1000);
N=1024;                                                       % ����Ƶ�ʺ�FFT�ĳ���
len=length(y);
time=(0:len-1)/fs;                                         % ʱ��̶�
figure(1), subplot 311; plot(time,y,'k');           % �����źŲ���
title('(a)�źŲ���'); axis([0 max(time) -1 1]);
ylabel('��ֵ'); xlabel(['ʱ��/s' 10]); grid;

nn=1:N/2; ff=(nn-1)*fs/N;                       % ����Ƶ�ʿ̶�
z=Nrceps(y);                                            %��ȡ����
figure(1), subplot 312; plot(time,z,'k');       % ��������ͼ
title('(b)�źŵ���ͼ'); axis([0 time(512) -0.2 0.2]); grid; 
ylabel('��ֵ'); xlabel(['��Ƶ��/s' 10]);

yc=cceps(y);
yn=icceps(yc);
figure(1), subplot 313; plot(time,yn,'k');     % ��������ͼ
title('(c)�ָ��ź�'); axis([0 max(time) -1 1]);
ylabel('��ֵ'); xlabel(['ʱ��/s' 10]); grid;     