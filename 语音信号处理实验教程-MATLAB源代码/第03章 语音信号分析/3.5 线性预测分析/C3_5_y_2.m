%ʵ��Ҫ�������̱Ƚ�LPCԤ��ϵ���ĸ�Ƶ����FFTƵ��
clear all; clc; close all;

[x,fs]=wavread('C3_5_y.wav');                    % ������������
L=240;                                                  	% ֡��
p=12;                                                       % LPC�Ľ���
y=x(8001:8000+L);                                   % ȡһ֡����
ar=lpc(y,p);                                                % ����Ԥ��任
nfft=512;                                                   % FFT�任����
W2=nfft/2;
m=1:W2+1;                                               % ��Ƶ�ʲ����±�ֵ
Y=fft(y,nfft);                                              % �����ź�y��FFTƵ��
Y1=lpcff(ar,W2);                                        % ����Ԥ��ϵ����Ƶ��
% ��ͼ
subplot 211; plot(y,'k');
title('(a)һ֡�����źŵĲ���'); ylabel('��ֵ'); xlabel('(a)')
subplot 212; 
plot(m,20*log10(abs(Y(m))),'k','linewidth',1.5); 
line(m,20*log10(abs(Y1)),'color',[.6 .6 .6],'linewidth',2)
axis([0 W2+1 -30 25]); ylabel('��ֵ/db');
legend('FFTƵ��','LPC��',3); xlabel(['����' 10 '(b)'])
title('(b)FFTƵ�׺�LPC�׵ıȽ�');
