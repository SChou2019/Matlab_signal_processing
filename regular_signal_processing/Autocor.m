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

xlabel('时间/s');ylabel('振幅');title('正弦信号时程曲线图');
subplot(2,2,2)
plot(a,b);

xlabel('时间/s');ylabel('振幅');title('正弦信号自相关函数');
subplot(2,2,3)
plot(t,x1);

xlabel('时间/s');ylabel('振幅');title('正弦信号加噪声时程曲线图');
subplot(2,2,4)
plot(d,c);

xlabel('时间/s');ylabel('振幅');title('正弦加噪声信号自相关函数');