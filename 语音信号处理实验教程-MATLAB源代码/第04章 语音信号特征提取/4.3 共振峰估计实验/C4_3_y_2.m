%LPC�ڲ巨�Ĺ�������
clear all; clc; close all;

fle='C4_3_y.wav';                            % ָ���ļ���
[x,fs]=wavread(fle);                        % ����һ֡�����ź� 
u=filter([1 -.99],1,x);                     % Ԥ����
wlen=length(u);                             % ֡��
p=12;                                       % LPC����
freq=(0:256)*fs/512;                        % Ƶ�ʿ̶�

[F,Bw,pp,U]=Formant_Interpolation(u,p,fs);          %LPC�ڲ巨�����
plot(freq,U,'k');
title('�������ݺ�������������');
xlabel('Ƶ��/Hz'); ylabel('��ֵ');
ll=length(F);                             % ��������
for k=1 : ll
    line([F(k) F(k)],[0 pp(k)],'color','k','linestyle','-.');    
end
legend('������','�����λ��')
fprintf('F =%5.2f   %5.2f   %5.2f   %5.2f\n',F)
fprintf('Bw=%5.2f   %5.2f   %5.2f   %5.2f\n',Bw)
