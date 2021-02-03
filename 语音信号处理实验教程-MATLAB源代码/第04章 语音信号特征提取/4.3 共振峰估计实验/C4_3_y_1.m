%ʵ��Ҫ��һ�����׷���������
clear all; clc; close all;

waveFile='C4_3_y.wav';               % �����ļ���
[x, fs, nbits]=wavread(waveFile);                 % ����һ֡����
u=filter([1 -.99],1,x);                                   % Ԥ����
wlen=length(u);                                          % ֡��
cepstL=6;                                                   % ��Ƶ���ϴ������Ŀ��
wlen2=wlen/2;               
freq=(0:wlen2-1)*fs/wlen;                          % ����Ƶ���Ƶ�ʿ̶�
u2=u.*hamming(wlen);		                      % �źżӴ�����
U=fft(u2);                                                 % ��ʽ(4-26)����
U_abs=log(abs(U(1:wlen2)));                     % ��ʽ(4-27)����
 [Val,Loc,spect]=Formant_Cepst(u2,cepstL);       % ����������Ƶ��
FRMNT=freq(Loc);                                 % ����������Ƶ��
subplot(211)
plot(freq,U_abs,'k'); 
xlabel('Ƶ��/Hz'); ylabel('��ֵ/dB');
title('(a)�źŶ�����X\_i(k)')
axis([0 4000 -6 2]); grid;
subplot(212)
plot(freq,spect,'k','linewidth',2); 
hold on
xlabel('Ƶ��/Hz'); ylabel('��ֵ/dB');
title('(b)�����ߺ͹����ֵ')
fprintf('%5.2f   %5.2f   %5.2f   %5.2f\n',FRMNT);
for k=1 : 4
    subplot(212)
    plot(freq(Loc(k)),Val(k),'kO','linewidth',2);
    line([freq(Loc(k)) freq(Loc(k))],[-6 Val(k)],'color','k',...
        'linestyle','-.','linewidth',2);
end
