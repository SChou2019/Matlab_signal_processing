%ʵ��Ҫ�����С��Ӳ��ֵ��������
clear all; clc; close all;

[xx, fs] = wavread('C5_4_y.wav');           % ���������ļ�
xx=xx-mean(xx);                         % ����ֱ������
x=xx/max(abs(xx));                      % ��ֵ��һ��
N=length(x);
%-------------------------����ָ��ǿ�ȵ�����---------------------------------
SNR=5;
s=awgn(x,SNR,'measured','db');               % ��������
wname='db7';

jN=6;  %�ֽ�Ĳ���
snrs=20*log10(norm(x)/norm(s-x));
signal=Wavelet_Hard(s,jN,wname);
signal=signal/max(abs(signal));
snr1=SNR_Calc(x,s);            % �����ʼ�����
snr2=SNR_Calc(x,signal);            % ���㽵���������
snr=snr2-snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n',snr1,snr2,snr);
% ��ͼ
time=(0:N-1)/fs;                        % ����ʱ��
subplot 311; plot(time,x,'k'); grid; axis tight;
title('����������'); ylabel('��ֵ')
subplot 312; plot(time,s,'k'); grid; axis tight;
title(['�������� �����=' num2str(SNR) 'dB']); ylabel('��ֵ')
subplot 313; plot(time,signal,'k');grid;%hold on;
title('�˲�����'); ylabel('��ֵ'); xlabel('ʱ��/s');
%--------------------------------------------------------------------------
