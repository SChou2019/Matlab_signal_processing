%  ʵ��Ҫ���������׷��������ڼ��
clc; close all; clear all;
wlen=320; inc=80;              % ��֡��֡����֡��
T1=0.05;                       % ���û����˵���Ĳ���
[x,fs]=wavread('C4_2_y.wav');                        % ����wav�ļ�
x=x-mean(x);                                % ��ȥֱ������
x=x/max(abs(x));                            % ��ֵ��һ��

[voiceseg,vosl,SF,Ef,period]=pitch_Ceps(x,wlen,inc,T1,fs); %���ڵ��׷��Ļ������ڼ��

fn=length(SF);
time = (0 : length(x)-1)/fs;                % ����ʱ������
frameTime = FrameTimeC(fn, wlen, inc, fs);  % �����֡��Ӧ��ʱ������
% ��ͼ
subplot 211, plot(time,x,'k');  title('�����ź�')
axis([0 max(time) -1 1]); ylabel('��ֵ');
subplot 212; plot(frameTime,period,'k');
xlim([0 max(time)]); title('���׷��������ڼ��'); 
xlabel('ʱ��/s'); ylabel('������');
for k=1 : vosl                              % ����л���
    nx1=voiceseg(k).begin;
    nx2=voiceseg(k).end;
    nxl=voiceseg(k).duration;
    fprintf('%4d   %4d   %4d   %4d\n',k,nx1,nx2,nxl);
    subplot 211
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
    subplot 212
    line([frameTime(nx1) frameTime(nx1)],[0 150],'color','r','linestyle','-');
    line([frameTime(nx2) frameTime(nx2)],[0 150],'color','b','linestyle','--');
end
