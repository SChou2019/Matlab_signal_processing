clear all;close all;clc
%%%%%%%%%%%%参数设置%%%%%%%%%%%%
s.M = 16; %阵元数
s.lamda = 1; %设置波长
s.d = 1/2*s.lamda; %取1/2波长为阵元距
s.SNR = 10; %信噪比
s.INR = 10; %干燥比
s.Noise = 1; %高斯噪声方差
s.SigDOA = [30 45]; %信号
s.IntDOA = [0 50 70]; %干扰
s.Theta = -90:1:90; %扫描范围
s.L = 200; %快拍数/采样数
s.fc = 100; %信号频率
s.fs=10000; %采样频率
%%%%%%%%%%生成信号%%%%%%%%%%%%
s.SigNum = length(s.SigDOA); %信号数
s.IntNum = length(s.IntDOA); %干扰数
s.ds = exp(1j*2*pi*(0:s.SigNum-1)'*s.fc/s.fs*(0:s.L-1)); %信号?
s.di = randn(s.IntNum,s.L)+1j*randn(s.IntNum,s.L); %干扰
s.As = exp(-1j*2*pi*s.d/s.lamda*(0:s.M-1)'*sin(s.SigDOA*pi/180)); %信号方向矢量
s.Ai = exp(-1j*2*pi*s.d/s.lamda*(0:s.M-1)'*sin(s.IntDOA*pi/180)); %干扰方向矢量
s.Xs = sqrt(10^(s.SNR/10))*s.As*s.ds; %构造信号
s.Xi = sqrt(10^(s.INR/10))*s.Ai*s.di; %构造干扰
s.Nos = (randn(s.M,s.L)+1j*randn(s.M,s.L))/sqrt(2); %构造噪声
s.Y = s.Xs+s.Xi+s.Nos; %合成阵列接收信号
SS = CBF(s);
figure(1)
plot(s.Theta,SS.P)
function s = CBF(s)
    s.X = s.Xs+s.Nos; %构造接收数据
    s.R=s.X*s.X'/s.L; %协方差处理
    a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta*pi/180)); %扫描
    s.P=a'*s.R*a; %常规波束形成输出功率
    s.P=diag(abs(s.P)); %提取特征值
    s.P=s.P/max(s.P); %输出功率归一化
    s.P=20*log10(abs(s.P));
end

function s=MMSE(s)
    % %%%%%%%%%%%数据处理%%%%%%%%%%%%
    s.X = s.Xs+s.Xi+s.Nos; %信号+干扰+噪声
    s.Rx=s.X*s.X'/s.L; %接收数据相关矩阵
    s.r_xd = s.X*s.ds'; %接收数据和期望信号互相关矩阵
    Wopt = s.Rx'*s.r_xd; %最佳权向量公式
    for t = 1:length(s.Theta)
    a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
    s.P(t) = Wopt'*a;
    end
    s.P=abs(s.P);
    s.P=s.P/max(max(s.P));
    s.P=20*log10(abs(s.P));
end