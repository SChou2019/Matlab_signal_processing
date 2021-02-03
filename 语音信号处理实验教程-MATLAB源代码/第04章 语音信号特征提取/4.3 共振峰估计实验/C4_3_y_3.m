% LPC������Ĺ�������
clear all; clc; close all;

fle='C4_3_y.wav';                            % ָ���ļ���
[xx,fs]=wavread(fle);                       % ����һ֡�����ź�
u=filter([1 -.99],1,xx);                    % Ԥ����
wlen=length(u);                             % ֡��
p=12;                                       % LPC����
n_frmnt=4;                                  % ȡ�ĸ������
freq=(0:256)*fs/512;                        % Ƶ�ʿ̶�
df=fs/512;                                  % Ƶ�ʷֱ���

[F,Bw,U]=Formant_Root(u,p,fs,n_frmnt);
plot(freq,U,'k');
title('�������ݺ�������������');
xlabel('Ƶ��/Hz'); ylabel('��ֵ/dB');
p1=length(F);                              % �ڹ���崦����
m=floor(F/df);
pp=U(m);                                    %��������
for k=1 : p1
    line([F(k) F(k)],[-5 pp(k)],'color','k','linestyle','-.');
end
legend('������','�����λ��')
fprintf('F0=%5.2f   %5.2f   %5.2f   %5.2f\n',F);
fprintf('Bw=%5.2f   %5.2f   %5.2f   %5.2f\n',Bw);

