%ʵ��Ҫ��һ���׼�����������
clear all; clc; close all;

[xx, fs] = wavread('C5_2_y.wav');           % ���������ļ�
xx=xx-mean(xx);                         % ����ֱ������
x=xx/max(abs(xx));                      % ��ֵ��һ��

IS=0.25;                                % ����ǰ���޻��γ���
wlen=200;                               % ����֡��Ϊ25ms
inc=80;                                 % ����֡��Ϊ10ms
SNR=5;                                  % ���������SNR
N=length(x);                            % �źų���
time=(0:N-1)/fs;                        % ����ʱ��

signal=awgn(x,SNR,'measured','db');               % ��������
snr1=SNR_Calc(x,signal);            % �����ʼ�����
NIS=fix((IS*fs-wlen)/inc +1);           % ��ǰ���޻���֡��

a=4; b=0.001;                           % ���ò���a��b
output=SpectralSub(signal,wlen,inc,NIS,a,b);% �׼�
snr2=SNR_Calc(x,output);            % �����׼���������
snr=snr2-snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n',snr1,snr2,snr);

% ��ͼ
subplot 311; plot(time,x,'k'); grid; axis tight;
title('����������'); ylabel('��ֵ')
subplot 312; plot(time,signal,'k'); grid; axis tight;
title(['�������� �����=' num2str(SNR) 'dB']); ylabel('��ֵ')
subplot 313; plot(time,output,'k');grid;%hold on;
title('�׼�����'); ylabel('��ֵ'); xlabel('ʱ��/s');



        
