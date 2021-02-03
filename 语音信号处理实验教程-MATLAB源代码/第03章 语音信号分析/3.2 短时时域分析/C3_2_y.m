%ʵ��Ҫ�󣺶�ʱʱ������������㲢��ʾ
clear all; clc; close all;

[x,Fs]=wavread('C3_2_y.wav');       % ���������ļ�

wlen=200; inc=100;          % ����֡����֡��
win=hanning(wlen);          % ����������
N=length(x);                    % �źų���
time=(0:N-1)/Fs;                % ������źŵ�ʱ��̶�
En=STEn(x,win,inc);             %��ʱ����
Mn=STMn(x,win,inc);             %��ʱƽ������
Zcr=STZcr(x,win,inc);              %��ʱ������
%�˴�������3��������ͬ�����صĲ����������Ǿ�����Ϊһ֡�źŵõ��Ĳ���һ����ֵ
X=enframe(x,win,inc)';     % ��֡
xn=X(:);
Ac=STAc(X);                         %�����ʱ�����
Ac=Ac(:);
Amdf=STAmdf(X);             %�����ʱ���Ȳ�
Amdf=Amdf(:);

fn=length(En);             % ���֡��

figure(1)
subplot 311; plot(time,x,'b'); axis tight% ����ʱ�䲨�� 
title('(a)��������');
ylabel('��ֵ'); xlabel(['ʱ��/s' 10 ]); 
frameTime=FrameTimeC(fn,wlen,inc,Fs);   % ���ÿ֡��Ӧ��ʱ��
subplot 312; plot(frameTime,Mn,'b')     % ������ʱ����ͼ
title('(b)��ʱ����');
ylabel('��ֵ'); xlabel(['ʱ��/s' 10 ]);  
subplot 313; plot(frameTime,En,'b')     % ������ʱ����ͼ
title('(c)��ʱ����');
 ylabel('��ֵ'); xlabel(['ʱ��/s' 10 '(b)']);

 figure(2)
subplot 211; plot(time,x,'b'); axis tight% ����ʱ�䲨�� 
title('(a)��������');
ylabel('��ֵ'); xlabel(['ʱ��/s' 10 ]); 
subplot 212; plot(frameTime,Zcr,'b')     % ������ʱ������ͼ
title('(b)��ʱ������');
ylabel('��ֵ'); xlabel(['ʱ��/s' 10 ]);  

figure(3)
subplot 211; plot(xn,'b'); % ����ʱ�䲨�� 
title('(a)��������');
ylabel('��ֵ'); xlabel(['����' 10 ]); 
subplot 212; plot(Ac,'b')                    % ������ʱ�����ͼ
title('(b)��ʱ�����');
ylabel('��ֵ'); xlabel(['����' 10 ]);  

figure(4)
subplot 211; plot(xn,'b'); % ����ʱ�䲨�� 
title('(a)��������');
ylabel('��ֵ'); xlabel(['����' 10 ]); 
subplot 212; plot(Amdf,'b')                     % ������ʱ���Ȳ�
title('(b)��ʱ���Ȳ�');
ylabel('��ֵ'); xlabel(['����' 10 ]);  
