%ʵ��Ҫ�����������ѹ������
%��д����ʵ��ʽ��2-4����Ҫ�������������Ǹ���ȼ�ʱ���ɵõ�����ȼ���Ӧ����ѹ�����ߣ���ʹ��plot����������ߵ���ʾ��
clc
clear all
phon=50;

[spl,freq]=iso226(phon);                    %������ѹ��

figure(1)
semilogx(freq,spl,':','color','k')
axis([20,20000,-10,130])
title('Phon=50')
xlabel('Ƶ��(Hz)')
ylabel('��ѹ����(dB)')
set(gca,'ytick',-10:10:130)
grid on
box off