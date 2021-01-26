%数据验证
clc
clear all
fs = 500;
w = 2*pi*200;
t = 0:1/fs:20;
%data_example = cos(fs*t+pi/6);
data_example = 2*cos(w*t+pi/6);
len = length(data_example);
nfft =1024;
overlap_rate = 0.67;
inrease_rate = 1 - overlap_rate;
m = floor((len-nfft)/(inrease_rate*nfft)+1);

%test
tmp_first = data_example(1:1024);
%乘上增长系数变成小数
%tmp_end = data_example(1+(m-1)*inrease_rate*nfft:(m-1)*inrease_rate*nfft+nfft);

A = [];
for i = 1:m
    start_index = floor(1+(m-1)*inrease_rate*nfft);
    end_index = floor((m-1)*inrease_rate*nfft+nfft);
    tem = data_example(start_index:end_index);
    fft_complex = fft(tem);
    fft_fix = abs(fft_complex) * 2 /length(tem);
    fft_value = fft_fix(1:length(tem)/2);
    A(i,:) =  fft_value;
end

duration = 1024/500;

%每个时间段以起点表示
increase_duration = inrease_rate *  duration ;
n = 1:27;
t = (n-1)*increase_duration;
delta = fs/nfft;
f = (1:length(tem)/2)*delta;

surf(f,t,A)





