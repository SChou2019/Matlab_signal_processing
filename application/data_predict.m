%���ݲ�����Ϊ25600����ֹƵ��12800
fs = 25600;
source = load('D:\software_learning\matlab_code\source.mat');
source_data = source.Audio_import_Front_outer2_6;
respond = load('D:\software_learning\matlab_code\respond.mat');
respond_data = respond.Audio_import_Front_outer2_7;

[num,txt,raw]=xlsread('D:\software_learning\matlab_code\NR_curve.xlsx');
fc = num(1:64,1)/12800;
mag = num(1:64,2);


%mag=[1 2 0.5 0.5 0.01 2.5 1 0.125 ]; %�����˲�������
N=32767;                  %�˲�������
b=fir2(N,fc,mag,50550,50);      %��ƺ������˲���

%freqz(b);              %����Ƶ����Ӧ����

length = length(source_data);
for i in 

