%ʵ��Ҫ�����LMS����
close all;clear all; clc; 

[s, fs, bits] = wavread('C5_1_y.wav');           % ���������ļ�
s=s-mean(s);                                % ����ֱ������
s=s/max(abs(s));                        % ��ֵ��һ
N=length(s);                                % ��������
time=(0:N-1)/fs;                        % ����ʱ��̶�
SNR=5;                                      % ���������
r1=awgn(s,SNR,'measured','db');
M=64;                                       % ���ãͺ�mu
mu=0.001;  
itr=length(r1);
snr1=SNR_Calc(s,r1);                    % �����ʼ�����
[y,W,e]=LMS(r1,s,M,mu,itr);
output=e/max(abs(e));                 	% LMS�˲����
snr2=SNR_Calc(s,output);            % �����˲���������
snr=snr2-snr1;
SN1=snr1; SN2=snr2; SN3=snr;
fprintf('snr1=%5.4f   snr2=%5.4f    snr=%5.4f\n',snr1,snr2,snr);
% ��ͼ
subplot 311; plot(time,s,'k'); ylabel('��ֵ') 
ylim([-1 1 ]); title('ԭʼ�����ź�');
subplot 312; plot(time,r1,'k'); ylabel('��ֵ') 
ylim([-1 1 ]); title('���������ź�');
subplot 313; plot(time,output,'k'); 
ylim([-1 1 ]); title('LMS�˲���������ź�');
xlabel('ʱ��/s'); ylabel('��ֵ')
