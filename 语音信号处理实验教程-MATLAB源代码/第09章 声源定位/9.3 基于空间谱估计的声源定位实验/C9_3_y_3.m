%ʵ��Ҫ��������ͬ�ռ��׹��Ƶ���Դ��λ�㷨�Ƚ�
clc;
clear all;
close all;
set(0,'defaultaxesfontsize',9);            %���������С
load s.mat;
wnd=256;
inc=128;

[Angle1]=Spectrum_Method('capon',s1,s2,wnd,inc,45);
subplot(311)
plot(Angle1,'*'),axis tight
title('Capon')
xlabel('֡��')
ylabel('�Ƕ�/��')

[Angle2]=Spectrum_Method('music',s1,s2,wnd,inc,45);
subplot(312)
plot(Angle2,'*'),axis tight
title('Music')
xlabel('֡��')
ylabel('�Ƕ�/��')

[Angle3]=Spectrum_Method('esprit',s1,s2,wnd,inc,45);
subplot(313)
plot(Angle3,'*'),axis tight
title('Esprit')
xlabel('֡��')
ylabel('�Ƕ�/��')