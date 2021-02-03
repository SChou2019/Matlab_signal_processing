%ʵ��Ҫ��һ����������������
clear all; clc; close all;
[x,fs,nbit]=wavread('C2_5_y_1.wav');     % ���������ļ�
len=length(x);
n=0.5:0.3/(len-1):0.8;                      %������������
x=x+n';                                         %��������������
t=(0:length(x)-1)/fs;                       % ����ʱ��
y=detrend(x);                               % ��������������
y=y/max(abs(y));                          % ��ֵ��һ��
subplot 211; plot(t,x,'k');               % ��������������������ź�x
title('��������������ź�');
xlabel('ʱ��/s'); ylabel('��ֵ');
subplot 212; plot(t,y,'k');               % ��������������������ź�y
xlabel('ʱ��/s'); ylabel('��ֵ');
title('����������������ź�');

%��������ʽ������
clear all; clc; 
[x,fs,nbit]=wavread('C2_5_y_1.wav');         % ����C2_5_y_1.wav�ļ�
len=length(x);
n=0:1/(len-1):1;
nn=n.^2-0.5;
x=x+nn';
[y,xtrend]=detrendN(x, fs, 2);          % ����detrendN����������
t=(0:length(x)-1)/fs;                        % ����ʱ��
figure
subplot 211; plot(t,x,'k');                 % ��������������������ź�x
line(t,xtrend,'color','r','linewidth',2); % ��������������
ylim([-1.5 1]);
title('��������������ź�');
legend('��������������ź�','�������ź�',4)
xlabel('ʱ��/s'); ylabel('��ֵ');
subplot 212; plot(t,y,'k');               % ��������������������ź�y
xlabel('ʱ��/s'); ylabel('��ֵ');
title('����������������ź�');


