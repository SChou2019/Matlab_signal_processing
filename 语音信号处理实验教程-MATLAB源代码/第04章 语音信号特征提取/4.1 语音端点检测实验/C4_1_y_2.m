%ʵ��Ҫ���������ط��˵���
clear all; clc; close all;
IS=0.25;                                % ����ǰ���޻��γ���
wlen=200;                               % ����֡��Ϊ25ms
inc=80;                                 % ��֡��

[xx,fs]=wavread('C4_1_y');                   % ��������
xx=xx-mean(xx);                         % ����ֱ������
x=xx/max(abs(xx));                      % ��ֵ��һ��
N=length(x);                            % ȡ�źų���
time=(0:N-1)/fs;                        % ����ʱ��
wnd=hamming(wlen);                      % ���ô�����
NIS=fix((IS*fs-wlen)/inc +1);           % ��ǰ���޻���֡��

% y=enframe(x,wnd,inc)';             % ��֡
% fn=size(y,2);                           % ��֡��
th1=1.1;
th2=1.3;
[voiceseg,vsl,SF,NF,Rum]=vad_corr(x,wnd,inc,NIS,th1,th2);% ����غ����Ķ˵���
fn=length(SF);
frameTime=FrameTimeC(fn, wlen, inc, fs);% �����֡��Ӧ��ʱ��
% ��ͼ
subplot 211; plot(time,x,'k');
title('����������');
ylabel('��ֵ'); axis([0 max(time) -1 1]);
subplot 212; plot(frameTime,Rum,'k');
title('��ʱ����غ���'); axis([0 max(time) 0 1]);
xlabel('ʱ��/s'); ylabel('��ֵ'); 
for k=1 : vsl                           % ��������˵�
    nx1=voiceseg(k).begin; nx2=voiceseg(k).end;
    subplot 211; 
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','LineStyle','--');
    subplot 212; 
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','LineStyle','--');
end



