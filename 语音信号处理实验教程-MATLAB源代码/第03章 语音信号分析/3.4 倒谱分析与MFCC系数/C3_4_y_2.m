%ʵ��Ҫ�����DCTϵ���ļ����벢�ָ�
clear all; clc; close all;

f=50;                                               % �ź�Ƶ��
fs=1000;                                        % ����Ƶ��
N=1000;                                     % ��������
n=0:N-1;
xn=cos(2*pi*f*n/fs);                    % ������������
y=dct(xn) ;                                     % ��ɢ���ұ任
num=find(abs(y)<5);                     % Ѱ�����ұ任���ֵС��5������
y(num)=0;                                    % �Է�ֵС��5������ķ�ֵ����Ϊ0
zn=idct(y);                                     % ��ɢ������任
subplot 211; plot(n,xn,'k');                % ����xn��ͼ
title('(a)ԭʼ�ź�'); xlabel(['����' 10 ]); ylabel('��ֵ');
subplot 212; plot(n,zn,'k');                % ����zn��ͼ
title('(b)�ؽ��ź�'); xlabel(['����' 10 ]); ylabel('��ֵ');
% �����ؽ���
rp=100-norm(xn-zn)/norm(xn)*100
