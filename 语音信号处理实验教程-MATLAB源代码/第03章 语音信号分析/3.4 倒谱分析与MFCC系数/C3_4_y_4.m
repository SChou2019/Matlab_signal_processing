%ʵ��Ҫ���ģ�MFCC�������
clear all; clc; close all;

[x1,fs]=wavread('C3_4_y_4.wav');                            % �����ź�C3_4_y_4.wav
wlen=200;                                                                  % ֡��
inc=80;                                                                     % ֡��
num=8;                                                                      %Mel�˲�������
x1=x1/max(abs(x1));                                                 % ��ֵ��һ��
time=(0:length(x1)-1)/fs;
subplot 211; plot(time,x1,'b') 
title('(a)�����ź�');
ylabel('��ֵ'); xlabel(['ʱ��/s' ]);  
ccc1=Nmfcc(x1,fs,num,wlen,inc);
fn=size(ccc1,1)+4;                                                  %ǰ�������֡������
cn=size(ccc1,2);
z=zeros(1,cn);
ccc2=[z;z;ccc1;z;z];
frameTime=FrameTimeC(fn,wlen,inc,fs);               % ���ÿ֡��Ӧ��ʱ��
subplot 212; plot(frameTime,ccc2(:,1:cn/2),'b')      % ����ÿͨ����MFCCϵ��
title('(b)MFCCϵ��');
ylabel('��ֵ'); xlabel(['ʱ��/s' ]);  

