% ʵ��Ҫ��һ���������ڼ��Ķ˵����㷨
clc; close all; clear all;
wlen=320; inc=80;              % ��֡��֡����֡��
T1=0.05;                       % ���û����˵���Ĳ���
[x,fs]=wavread('C4_2_y.wav');                        % ����wav�ļ�
x=x-mean(x);                                % ��ȥֱ������
x=x/max(abs(x));                            % ��ֵ��һ��

[voiceseg,vosl,SF,Ef]=pitch_vad(x,wlen,inc,T1);   % �����Ķ˵���
fn=length(SF);
time = (0 : length(x)-1)/fs;                % ����ʱ������
frameTime = FrameTimeC(fn, wlen, inc, fs);  % �����֡��Ӧ��ʱ������
plot(time,x,'k');  title('�����ź�')
axis([0 max(time) -1 1]); ylabel('��ֵ');
xlabel('ʱ��/s');
for k=1 : vosl                              % ����л���
    nx1=voiceseg(k).begin;
    nx2=voiceseg(k).end;
    nxl=voiceseg(k).duration;
    fprintf('%4d   %4d   %4d   %4d\n',k,nx1,nx2,nxl);
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
end
