%ʵ��Ҫ���������Դ˫��˷�ķ���弤��Ӧʵ��
clear all;
close all;
fs=8000;                                                                    %����Ƶ��
path_length = 256;                                                      %·������
mic1=[0.9 1.5 1.5];                                                       %��˷�1λ��
mic2=[1.7 1.5 1.5];                                                           %��˷�2λ��
n=12;                                                                          %����Դ����
r=0.25;                                                                          %����ϵ��
c=340;                                                                          %����
rm=[4 4 3];                                                                   %����ߴ�
src=[2.1 2.5 1.5];                                                           %��Դλ��
h=rir(fs, mic1, n, r, rm, src);
h1=h(1:path_length);

h=rir(fs, mic2, n, r, rm, src);
h2=h(1:path_length);
figure(1);
subplot(2,1,1), plot(h1),axis([1,path_length,-1,1]);
ylabel('����')
xlabel('����')
title('��˷�1�ĳ弤��Ӧ')
subplot(2,1,2), plot(h2),axis([1,path_length,-1,1]);
ylabel('����')
xlabel('����')
title('��˷�2�ĳ弤��Ӧ')
save ('h.mat', 'h1', 'h2');

