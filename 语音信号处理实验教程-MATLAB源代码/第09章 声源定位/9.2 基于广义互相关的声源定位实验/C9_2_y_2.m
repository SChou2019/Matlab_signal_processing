%ʵ��Ҫ��ʱ�ӹ����㷨�Ա�ʵ��
clear all;
close all;
clc
load s.mat;
s1=s1_10db(:)';             %��Ϊ������
s2=s2_10db(:)';             %��Ϊ������
wnd=512;
inc=256;

[delay]=GCC_Method('standard',s1,s2,wnd,inc);
subplot(411)
plot(delay-wnd,'*')
ylim([0,12])
title('��׼GCC')
xlabel('֡��')
ylabel('��ʱ/��')

[delay]=GCC_Method('phat',s1,s2,wnd,inc);
subplot(412)
plot(delay-wnd,'*')
ylim([0,12])
title('Phat-GCC')
xlabel('֡��')
ylabel('��ʱ/��')

[delay]=GCC_Method('scot',s1,s2,wnd,inc);
subplot(413)
plot(delay-wnd,'*')
ylim([0,12])
title('Scot-GCC')
xlabel('֡��')
ylabel('��ʱ/��')

[delay]=GCC_Method('ml',s1,s2,wnd,inc);
 subplot(414)
plot(delay-wnd,'*')
ylim([0,12])
title('Ml-GCC')
xlabel('֡��')
ylabel('��ʱ/��')


    
