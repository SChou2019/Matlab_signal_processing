%ʵ��Ҫ��һ��˫���޷��˵���
clear all; clc; close all;

[x,fs]=wavread('C4_1_y.wav');                    % ���������ļ�
x=x/max(abs(x));                        % ���ȹ�һ��
N=length(x);                            % ȡ�źų���
time=(0:N-1)/fs;                        % ����ʱ��
subplot 311
plot(time,x,'k');         
title('˫���޷��Ķ˵���');
ylabel('��ֵ'); axis([0 max(time) -1 1]); 
xlabel('ʱ��/s');
wlen=200; inc=80;                       % ��֡����
IS=0.1; overlap=wlen-inc;               % ����IS
NIS=fix((IS*fs-wlen)/inc +1);           % ����NIS
fn=fix((N-wlen)/inc)+1;                 % ��֡��
frameTime=FrameTimeC(fn, wlen, inc, fs);% ����ÿ֡��Ӧ��ʱ��
[voiceseg,vsl,SF,NF,amp,zcr]=vad_TwoThr(x,wlen,inc,NIS);  % �˵���
subplot 312
plot(frameTime,amp,'k');         
ylim([min(amp) max(amp)])
title('��ʱ����');
ylabel('��ֵ'); 
xlabel('ʱ��/s');
subplot 313
plot(frameTime,zcr,'k');     
ylim([min(zcr) max(zcr)])
title('��ʱ������');
ylabel('��ֵ'); 
xlabel('ʱ��/s');
for k=1 : vsl                           % ������ֹ��λ��
    subplot 311
    nx1=voiceseg(k).begin; nx2=voiceseg(k).end;
    nxl=voiceseg(k).duration;
    line([frameTime(nx1) frameTime(nx1)],[-1.5 1.5],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[-1.5 1.5],'color','b','LineStyle','--');
    subplot 312
    line([frameTime(nx1) frameTime(nx1)],[min(amp) max(amp)],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[min(amp) max(amp)],'color','b','LineStyle','--');    
    subplot 313
    line([frameTime(nx1) frameTime(nx1)],[min(zcr) max(zcr)],'color','r','LineStyle','-');
    line([frameTime(nx2) frameTime(nx2)],[min(zcr) max(zcr)],'color','b','LineStyle','--');    
end
