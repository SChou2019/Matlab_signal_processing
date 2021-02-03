clc
clear all
%***************1.���Ҳ�****************%
fs =100;                                            %�趨����Ƶ��
N =128;
n =0:N -1;
t = n/ fs;
f0 =10;                                             %�趨�����ź�Ƶ��
%���������ź�
x = sin(2*pi*f0*t);
figure(1);
subplot(231);
plot(t,x); %�������źŵ�ʱ����
xlabel('ʱ��/ s');
ylabel('��ֵ');
title('ʱ����');
grid;
%����FFT �任����Ƶ��ͼ
y = fft(x,N);                                       %����FFT �任
mag = abs(y);                                   %���ֵ
f = (0:length(y) -1)'*fs/ length(y);        %���ж�Ӧ��Ƶ��ת��
subplot(232);
plot(f,mag);                                        %��Ƶ��ͼ
axis([0,100,0,80]);
xlabel('Ƶ��/ Hz');
ylabel('��ֵ');
title('��Ƶ��ͼ');
grid;
%���������
sq = abs(y);
subplot(233);
plot(f,sq);
xlabel('Ƶ��/ Hz');
ylabel('��������');
title('��������');
grid;
%������
power = sq.^2;
subplot(234);
plot(f,power);
xlabel('Ƶ��/ Hz');
ylabel('������');
title('������');
grid;
%�������
ln = log(sq);
subplot(235);
plot(f,ln);
xlabel('Ƶ��/ Hz');
ylabel('������');
title('������');
grid;
%��IFFT �ָ�ԭʼ�ź�
xifft = ifft(y);
magx = real(xifft);
ti = [0:length(xifft)-1] / fs;
subplot(236);
plot(ti,magx);
xlabel('ʱ��/ s');
ylabel('��ֵ');
title('IFFT ����źŲ���');
grid;
%****************2.������****************%
fs =50;                                              %�趨����Ƶ��
t = -5:0.1:5;
x = rand(1,100);
figure(2);
subplot(231);
plot(t(1:100),x);                               %����������ʱ����
xlabel('ʱ��(s)');
ylabel('��ֵ');
title('ʱ����');
grid;
%����FFT �任����Ƶ��ͼ
y = fft(x);                                     %����FFT �任
mag = abs(y);                               %���ֵ
f = (0:length(y) -1)'*fs/ length(y); %���ж�Ӧ��Ƶ��ת��
subplot(232);
plot(f,mag);                                    %��Ƶ��ͼ
xlabel('Ƶ��/ Hz');
ylabel('��ֵ');
title('��Ƶ��ͼ');
grid;
%���������
sq = abs(y);
subplot(233);
plot(f,sq);
xlabel('Ƶ��/ Hz');
ylabel('��������');
title('��������');
grid;
%������
power = sq.^2;
subplot(234);
plot(f,power);
xlabel('Ƶ��/ Hz');
ylabel('������');
title('������');
grid;
%�������
ln = log(sq);
subplot(235);
plot(f,ln);
xlabel('Ƶ��/ Hz');
ylabel('������');
title('������');
grid;
%��IFFT �ָ�ԭʼ�ź�
xifft = ifft(y);
magx = real(xifft);
ti = [0:length(xifft)-1] / fs;
subplot(236);
plot(ti,magx);
xlabel('ʱ��/ s');
ylabel('��ֵ');
title('IFFT ����źŲ���');
grid;