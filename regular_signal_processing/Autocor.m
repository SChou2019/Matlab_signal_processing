clear;
clc;
fs = 1000;
N = 1000;
n = 0:N-1;
t = n/fs;
x=sin(2*pi*15*t);

x1=awgn(x,0);

figure

subplot

[b,a]=xcorr(x,'unbiased');

[c,d]=xcorr(x1,'unbiased');
subplot(2,2,1)
plot(t,x);

xlabel('ʱ��/s');ylabel('���');title('�����ź�ʱ������ͼ');
subplot(2,2,2)
plot(a,b);

xlabel('ʱ��/s');ylabel('���');title('�����ź�����غ���');
subplot(2,2,3)
plot(t,x1);

xlabel('ʱ��/s');ylabel('���');title('�����źż�����ʱ������ͼ');
subplot(2,2,4)
plot(d,c);

xlabel('ʱ��/s');ylabel('���');title('���Ҽ������ź�����غ���');