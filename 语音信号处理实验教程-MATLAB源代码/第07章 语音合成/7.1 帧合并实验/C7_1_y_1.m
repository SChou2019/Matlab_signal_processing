%ʵ��Ҫ��һ�������ص���ӷ����źŷ�֡�뻹ԭ
clc
clear all
close all
[s,fs]=wavread('C7_1_y.wav');
winlen=256;
win=hamming(winlen);
overlap=100;
f=enframe(s,win,overlap);
fn=Filpframe_OverlapA(f,win,overlap);
subplot(211)
plot(s/max(abs(s)))
xlabel('����')
ylabel('����')
title('(a)ԭʼ�ź�')
subplot(212)
plot(fn/max(abs(fn)))
xlabel('����')
ylabel('����')
title('(b)��ԭ�ź�')