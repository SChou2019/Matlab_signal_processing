%ʵ��Ҫ��һ������Ԥ��ϵ���Ա�
clear all; clc; close all;

[x,fs]=wavread('C3_5_y.wav');    % ������������
L=240;                                      % ֡��
y=x(8001:8000+L);                   % ȡһ֡����
p=12;                                       % LPC�Ľ���
ar1=lpc(y,p);                            % MATLAB�Դ�������������Ԥ��任
ar2=lpc_coeff(y,p);                  % ��д�ĺ�����������Ԥ��任
est_x1=filter([0 -ar1(2:end)],1,y);       % ��LPC��Ԥ�����ֵ
est_x2=filter([0 -ar2(2:end)],1,y);       % �ñ�д������Ԥ�����ֵ
err1=y-est_x1;                            % LPC��Ԥ�����
err2=y-est_x2;                            % ��д������Ԥ�����

subplot 321; plot(x,'k'); axis tight;
title('(a)Ԫ��/a/����'); ylabel('��ֵ')
subplot 322; plot(y,'k'); xlim([0 L]); 
title('(b)һ֡����'); ylabel('��ֵ')
subplot 323; plot(est_x1,'k'); xlim([0 L]); 
title('(c)LPCԤ��ֵ'); ylabel('��ֵ')
subplot 324; plot(est_x2,'k'); xlim([0 L]); 
title('(d)lpc\_coeffԤ��ֵ'); ylabel('��ֵ')
subplot 325; plot(err1,'k'); xlim([0 L]); 
title('(e)LPCԤ�����'); ylabel('��ֵ'); xlabel('����')
subplot 326; plot(err2,'k'); xlim([0 L]); 
title('(f)lpc\_coeffԤ�����'); ylabel('��ֵ'); xlabel('����')









