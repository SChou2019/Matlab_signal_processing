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
overlap=wlen-inc;                       % ���ص�������
NIS=fix((IS*fs-wlen)/inc +1);           % ��ǰ���޻���֡��

th1=0.99;
th2=0.96;
[voiceseg,vsl,SF,NF,Enm]=vad_specEn(x,wnd,inc,NIS,th1,th2,fs);  % ���ط�  
fn=length(SF);
frameTime=FrameTimeC(fn, wlen, inc, fs);% �����֡��Ӧ��ʱ��
subplot 211; 
plot(time,x,'k'); hold on
title('��������');
ylabel('��ֵ'); axis([0 max(time) -1 1]);
subplot 212; plot(frameTime,Enm,'k');  
ylim([min(Enm) max(Enm)])
title('��ʱ�Ľ��Ӵ�����'); xlabel('ʱ��/s'); ylabel('����ֵ');
for k=1 : vsl                         
    nx1=voiceseg(k).begin; nx2=voiceseg(k).end;
    fprintf('%4d   %4d   %4d\n',k,nx1,nx2);
    subplot 211
    line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','LineStyle','--');
    subplot 212
    line([frameTime(nx1) frameTime(nx1)],[min(Enm) max(Enm)],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[min(Enm) max(Enm)],'color','b','LineStyle','--');
end

