%ʵ��Ҫ�����������źŲ���Ƶ�ʱ任
clc
clear all
[x,fs1]=wavread('C2_2_y.wav');
s1=1:length(x);
t1=s1/fs1;
xmax=max(abs(x));
x=x/xmax;

figure(1)
subplot(311)
plot(t1,x);
xlabel('ʱ��/s');
ylabel('��һ����ֵ');
title('(a)ԭʼ�ź�');

p=2;q=1;
x1=resample(x,p,q);
x1max=max(abs(x1));
x1=x1/x1max;
fa=fs1*p/q;
ta=(1:length(x1))/fa;
subplot(312);
plot(ta,x1);
xlabel('ʱ��/s');
ylabel('��һ����ֵ');
title('(b)2��������');

p=1;q=2;
x2=resample(x,p,q);
x2max=max(abs(x2));
x2=x2/x2max;
fb=fs1*p/q;
tb=(1:length(x2))/fb;
subplot(313);
plot(tb,x2);
xlabel('ʱ��/s');
ylabel('��һ����ֵ');
title('(c)1/2������');