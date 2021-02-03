%ʵ��Ҫ�󣺻�����������Ϣ����ʵ��
clc;
clear all;
close all;
[x,fs]=wavread('C8_2_y.wav');
message = zeros(1,100);
for i=1:100
    if (rand(1)>0.5)
        message(1,i)=1;
    else
        message(1,i)=0;
    end
end
N = 0.1*fs;
m0 = 0.4*N;
m1=0.2*N;
x_embeded=hide_EchoEmbed(x,message,N,m0,m1,0.7);
len  = min(length(message),length(x)/N);
mess = hide_EchoExtract( x_embeded,N,m0,m1,len);

figure;plot(x);title('ԭʼ�����ź�');
figure;plot(x_embeded);title('��Ƕ����Ϣ�������ź�');
figure;plot(mess,'*');title('��ȡ��������Ϣ');grid on;