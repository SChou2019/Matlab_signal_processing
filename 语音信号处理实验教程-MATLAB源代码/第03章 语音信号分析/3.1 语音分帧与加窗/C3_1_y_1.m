%ʵ��Ҫ��һ��������֡��ʾ
clc
clear all
close all

[x,Fs]=wavread('C3_1_y.wav');       % ���������ļ�
wlen=200; inc=100;          % ����֡����֡��
N=length(x);                    % �źų���
time=(0:N-1)/Fs;                % ������źŵ�ʱ��̶�
signal=enframe(x,wlen,inc)';     % ��֡

i=input('��������ʼ֡��(i):');
tlabel=i;
subplot 411; plot((tlabel-1)*inc+1:(tlabel-1)*inc+wlen,signal(:,tlabel),'b'); axis tight% ����ʱ�䲨�� 
xlim([(i-1)*inc+1 (i+2)*inc+wlen])
ylim([-0.1,0.1])
title(['(a)��ǰ����֡�ţ�', num2str(i)]);
ylabel('��ֵ'); xlabel('֡��'); 
tlabel=i+1;
subplot 412; plot((tlabel-1)*inc+1:(tlabel-1)*inc+wlen,signal(:,tlabel),'b'); axis tight% ����ʱ�䲨�� 
xlim([(i-1)*inc+1 (i+2)*inc+wlen])
ylim([-0.1,0.1])
title(['(b)��ǰ����֡�ţ�', num2str(i+1)]);
ylabel('��ֵ'); xlabel('֡��'); 
tlabel=i+2;
subplot 413; plot((tlabel-1)*inc+1:(tlabel-1)*inc+wlen,signal(:,tlabel),'b'); axis tight% ����ʱ�䲨�� 
xlim([(i-1)*inc+1 (i+2)*inc+wlen])
ylim([-0.1,0.1])
title(['(c)��ǰ����֡�ţ�', num2str(i+2)]);
ylabel('��ֵ'); xlabel('֡��'); 
tlabel=i+3;
subplot 414; plot((tlabel-1)*inc+1:(tlabel-1)*inc+wlen,signal(:,tlabel),'b'); axis tight% ����ʱ�䲨�� 
xlim([(i-1)*inc+1 (i+2)*inc+wlen])
ylim([-0.1,0.1])
title(['(d)��ǰ����֡�ţ�', num2str(i+3)]);
ylabel('��ֵ'); xlabel('֡��'); 
