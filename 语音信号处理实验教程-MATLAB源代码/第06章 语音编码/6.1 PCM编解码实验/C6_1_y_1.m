%ʵ��Ҫ��һ��PCM�����ʵ��
clear all;
close all;
[x fs numbits]= wavread('C6_1_y.wav');  
v=1;
xx=x/v;
sxx=floor(xx*4096);
y=pcm_encode(sxx);
yy=pcm_decode(y,v)';

nq=sum((x-yy).*(x-yy))/length(x);
sq=mean(yy.^2);
snr=(sq/nq);
t=(1:length(x))/fs;
subplot(211)
plot(t,x/max(abs(x)))
axis tight
title('(a)����ǰ����')
xlabel('ʱ��/s')
ylabel('����')
subplot(212)
plot(t,yy/max(abs(yy)))
axis tight
title('(b)���������')
xlabel('ʱ��/s')
ylabel('����')
snrq=10*log10(mean(snr))


