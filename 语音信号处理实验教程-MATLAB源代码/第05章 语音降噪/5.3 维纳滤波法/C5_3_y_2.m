%ʵ��Ҫ�����������������ȵ�ά���˲��㷨��������
clear all; clc; close all;

[xx, fs] = wavread('C5_3_y.wav');           % ���������ļ�
xx=xx-mean(xx);                         % ����ֱ������
x=xx/max(abs(xx));                      % ��ֵ��һ��
IS=0.25;                                % ����ǰ���޻��γ���
wlen=200;                               % ����֡��Ϊ25ms
inc=80;                                 % ����֡��Ϊ10ms
SNR=5;                                  % ���������SNR
NIS=fix((IS*fs-wlen)/inc +1);           % ��ǰ���޻���֡��
alpha=0.95;

signal=awgn(x,SNR,'measured','db');               % ��������
output=Weina_Im(x,wlen,inc,NIS,alpha) ;
output=output/max(abs(output));
len=min(length(output),length(x));
x=x(1:len);
signal=signal(1:len);
output=output(1:len);

snr1=SNR_Calc(x,signal);            % �����ʼ�����
snr2=SNR_Calc(x,output);            % ���㽵���������
snr=snr2-snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n',snr1,snr2,snr);

% ��ͼ
time=(0:len-1)/fs;                        % ����ʱ��
subplot 311; plot(time,x,'k'); grid; axis tight;
title('����������'); ylabel('��ֵ')
subplot 312; plot(time,signal,'k'); grid; axis tight;
title(['�������� �����=' num2str(SNR) 'dB']); ylabel('��ֵ')
subplot 313; plot(time,output,'k');grid;%hold on;
title('�˲�����'); ylabel('��ֵ'); xlabel('ʱ��/s');



        
