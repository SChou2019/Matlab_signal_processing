%ʵ��Ҫ��������̱Ƚ�LPCC
clear all; clc; close all;

[x,fs]=wavread('C3_5_y.wav');            % ������������
L=240;                                              % ֡��
p=12;                                               % LPC�Ľ���
y=x(8001:8000+L);                           % ȡһ֡����
ar=lpc(y,p);                                       % ����Ԥ��任
lpcc1=lpc_lpccm(ar,p,p);
lpcc2=rceps(ar); 

subplot 211; plot(lpcc1(1:2:end),'k');
title('(a)����Ԥ��ϵ����LPCC'); ylabel('��ֵ'); xlabel(['����' ])
subplot 212; plot(lpcc2(1:p/2),'k');
title('(b)ֱ����LPCC'); ylabel('��ֵ'); xlabel(['����' ])
