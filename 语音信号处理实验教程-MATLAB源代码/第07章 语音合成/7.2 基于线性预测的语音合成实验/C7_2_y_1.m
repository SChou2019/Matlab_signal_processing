%ʵ��Ҫ��һ����������Ԥ��ϵ����Ԥ�����������ϳ�ʵ��
clear all; clc; close all;

[x, fs, bits] = wavread('C7_2_y.wav');           % ���������ļ�
x=x-mean(x);                            % ����ֱ������
x=x/max(abs(x));                        % ��ֵ��һ
xl=length(x);                           % ���ݳ���
time=(0:xl-1)/fs;                       % �����ʱ��̶�
p=12;                                   % LPC�Ľ���Ϊ12
wlen=200; inc=80;                       % ֡����֡��
msoverlap = wlen - inc;                 % ÿ֡�ص����ֵĳ���
y=enframe(x,wlen,inc)';                 % ��֡
fn=size(y,2);                           % ȡ֡��
% ��������:��ÿһ֡��LPCϵ����Ԥ�����
for i=1 : fn                            
    u=y(:,i);                           % ȡ��һ֡
    A=lpc(u,p);                         % LPC���ϵ��
    aCoeff(:,i)=A;                      % �����aCoeff������
    errSig = filter(A,1,u);             % ����Ԥ���������
    resid(:,i) = errSig;                % �����resid������
end
% �����ϳ�:��ÿһ֡�ĺϳ��������ӳ����������ź�
for i=1:fn                              
    A = aCoeff(:,i);                    % ȡ�ø�֡��Ԥ��ϵ��
    residFrame = resid(:,i);            % ȡ�ø�֡��Ԥ�����
    synFrame(i,:) = filter(1, A', residFrame); % Ԥ������,�ϳ�����
end;
outspeech=Filpframe_OverlapS(synFrame,wlen,inc);
ol=length(outspeech);
if ol<xl                                % ��outspeech����,ʹ��x�ȳ�
    outspeech=[outspeech zeros(1,xl-ol)];
else
    outspeech=outspeech(1:xl);
end

% ��ͼ
subplot 211; plot(time,x,'k');
xlabel(['ʱ��/s']); ylabel('��ֵ'); ylim([-1 1.1]);
title('(a)ԭʼ�����ź�')
subplot 212; plot(time,outspeech,'k');
xlabel(['ʱ��/s' ]); ylabel('��ֵ'); ylim([-1 1.1]);
title('(b)�ϳɵ������ź�')


