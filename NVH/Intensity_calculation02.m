clc;clear all;

channel_1 = load("C:\Users\Administrator\Desktop\export\channel01.mat");%ͨ��1��Time block����
channel_2 = load("C:\Users\Administrator\Desktop\export\channel02.mat");
channel_1_01 = channel_1.ExampleIntensity_Mainframe_(:,1);%ͨ��1��һ���������ݿ�
channel_2_01 = channel_2.ExampleIntensity_Mainframe_(:,1);%ͨ��1��һ���������ݿ�
Fs = 16384;%����Ƶ��
delta_r = 0.012;%̽ͷ���������
air_density = 1.239;%�����ܶ�

%Ƶ�� Intensity calculation

FFT_1 = fft(channel_1_01);%ͨ��1Ƶ�׼���
FFT_1 = FFT_1*2/length(channel_1_01);%��ֵ����
FFT_2 = fft(channel_2_01);
FFT_2 = FFT_2*2/length(channel_2_01);

delatf =  Fs/length(channel_1_01);%Ƶ�ʷֱ���
CPS = conj(FFT_1).*FFT_2;%����

CPS_f = CPS(1:floor(length(channel_1_01)/2)); 
CPS_imag = imag(CPS_f);%ȡ�鲿
f = (1:floor(length(channel_1_01)/2))*delatf;  %Ƶ�ʵ�
w_f = 2*pi*f';%ԲƵ��
denominator = air_density*w_f*delta_r;
%Intensity = -(imag(CPS)./(air_density*w_f*delta_r));
Intensity = -CPS_imag./denominator;
%plot��ʾ�������Ĵ���
figure(1)
plot(f,Intensity)

% octave plot
% ����Ƶ�����
Fc = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 ...
500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 ... 
6300 8000 10000 12500 16000 20000];% frequenze centrali ANSI
Fc = Fc(Fc>100);%��125Hz����Ƶ��
Fc = Fc(Fc<Fs/2);%���޵���Fs/2��Ҳ����ǿ������߶����
fl = Fc*2^(-1/6); % lower freq
fu = Fc*2^(1/6); % upper freq
fu(fu>Fs/2) = Fs/2 -1; 
numBands = length(Fc);
y = zeros(numBands,1);
for i = 1:numBands
    for j = 1:length(f)
        if f(j) >= fl(i) && f(j) <= fu(i)
            %ֱ�����
            y(i) =  y(i) + Intensity(j);
        end
    end
end

%�����Ĵ���
y_dB = 10*log10(abs(y)/1e-12);
%y_dB = abs(10*log10(y/1e-12));
figure(2)
bar(Fc,y_dB)