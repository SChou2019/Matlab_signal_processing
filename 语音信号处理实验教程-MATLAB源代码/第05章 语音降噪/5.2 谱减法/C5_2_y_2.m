%ʵ��Ҫ�����Boll�ĸĽ��׼���
clear all; clc; close all;

[xx,fs]=wavread('C5_2_y.wav');                   % ���������ļ�
xx=xx-mean(xx);                         % ����ֱ������
x=xx/max(abs(xx));                      % ��ֵ��һ��
SNR=10;                                 % ���������
signal=awgn(x,SNR,'measured','db');                % ��������
snr1=SNR_Calc(x,signal);            % �������������������
N=length(x);                            % �źų���
time=(0:N-1)/fs;                        % ����ʱ��̶�
IS=.15;                                 % % ����ǰ���޻��γ���
wlen=200;                               % ����֡��Ϊ25ms
inc=80;                                 % ����֡��Ϊ10ms
Gamma=1;                                %���ȼ�Ȩ���Ľ��׼����Ĳ�����
Beta=.03;
NIS=fix((IS*fs-wlen)/inc +1);           % ��ǰ���޻���֡��
output=SpectralSubIm(signal,wlen,inc,NIS,Gamma,Beta);          % ����SSBoll79�������׼�
output=output/max(abs(output));
ol=length(output);                      % ��output������x�ȳ�
if ol<N
    output=[output; zeros(N-ol,1)];
end
snr2=SNR_Calc(x,output);            % �����׼���������
snr=snr2-snr1;
fprintf('snr1=%5.4f   snr2=%5.4f   snr=%5.4f\n',snr1,snr2,snr);

% ��ͼ
subplot 311; plot(time,x,'k'); grid; axis tight;
title('����������'); ylabel('��ֵ')
subplot 312; plot(time,signal,'k'); grid; axis tight;
title(['�������� �����=' num2str(SNR) 'dB']); ylabel('��ֵ')
subplot 313; plot(time,output,'k');grid; ylim([-1 1]);
title('�׼�����'); ylabel('��ֵ'); xlabel('ʱ��/s');