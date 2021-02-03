%ʵ��Ҫ��һ��LPC�����ʵ��
clc
close all
clear all
[x,fs,nbits]=wavread('C6_2_y');     
LENGTH=length(x);                          
n=0:1/fs:(LENGTH-1)/fs;   

% ---------------------------------------------------------------------- %
% ԭʼ�����ź� %
subplot(2,1,1),plot(n*1000, x),grid ,hold on
xlabel('ʱ��/��');    ylabel('����');

% �򵥱������ʾ %
order=12;                                               % ����
[a,g]=lpc(x,order);                                   % Ԥ��ϵ�������������
est_x=filter([0 -a(2:end)],1,x);                  % ������ź�
plot(n*1000,est_x,'r--'),hold off
title('ԭʼ�ź��������ź�')
legend('ԭʼ�ź�','������ź�')

% ---------------------------------------------------------------------- %
error=x-est_x;                                    % Ԥ�������ʾ
subplot(2,1,2),plot(n*1000,error), grid;
xlabel('ʱ��/��');    ylabel('����');
title('Ԥ�����')