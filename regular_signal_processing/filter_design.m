%���ݲ�����Ϊ25600����ֹƵ��10000
source_data = load('D:\software_learning\matlab_code\source.mat');
respond_data = load('D:\software_learning\matlab_code\respond.mat');
%[num,txt,raw]=xlsread('D:\software_learning\matlab_code\test.xlsx');
[num,txt,raw]=xlsread('D:\software_learning\matlab_code\test.xlsx');
fc = num(6:1606,1)/1600;
mag = num(6:1606,2);
%fc=[0 0.125 0.25 0.35 0.45 0.55 0.75 1];   %��ֹƵ��

%mag=[1 2 0.5 0.5 0.01 2.5 1 0.125 ]; %�����˲�������
N=100;                  %�˲�������
b=fir2(N,fc,mag,51);      %��ƺ������˲���

freqz(b);              %����Ƶ����Ӧ����

ref = 2/(power(10,6/20));% should be 1
mag_db = 20*log10(mag/ref);