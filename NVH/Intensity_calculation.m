clc;clear all;
data = load("C:\Users\Administrator\Desktop\export\Throughput.mat");
channel_1 = load("C:\Users\Administrator\Desktop\export\channel01.mat");
channel_2 = load("C:\Users\Administrator\Desktop\export\channel02.mat");
channel_1_01 = channel_1.ExampleIntensity_Mainframe_(:,1);%ͨ��1��һ���������ݿ�
channel_2_01 = channel_2.ExampleIntensity_Mainframe_(:,1);%ͨ��1��һ���������ݿ�
fs = 16384;
delta_r = 0.012;%���������
air_density = 1.239;
% source = data.BSWA_2s___1_1(:,1);
% respond = data.BSWA_2s___1_1(:,2);

%Intensity calculation
%��ֵ��У������
FFT_s = fft(source);
FFT_s = FFT_s*2/length(source);
FFT_r = fft(respond);
FFT_r = FFT_r*2/length(source);

CPS = conj(FFT_s).*FFT_r;
CPS_f = CPS(1:floor(length(source)/2)); 
CPS_f_imag = imag(CPS_f );
f = (1:floor(length(source)/2))*fs/length(source);  %Ƶ�ʵ�
w_f = 2*pi*f';
denominator = air_density*w_f*delta_r;
%�����
%Intensity = -(imag(CPS)./(air_density*w_f*delta_r));
Intensity = -CPS_f_imag./denominator;