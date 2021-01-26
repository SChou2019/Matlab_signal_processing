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

%常规波束形成是对不含干扰的信号进行波束形成， 构造接收数据就是构造信号加上噪声；
%常规波束形成适用于单个信号且不含干扰的声源信号波束形成。 瑞利准则说
%明常规波束形成法的故有缺点就是角分辨率低， 要提高角分变率就需要增大阵元
%间隔或增加阵元个数。
function s = CBF(s)
s.X = s.Xs + s.Nos;                                       %构造接收数据
s.R = s.X*s.X'/s.L;                                       %协方差处理
a = exp(-1i*pi*(0:s.M-1)*sin(s.Theta*pi/180));            %扫描
s.P = a'*s.R*a;                                           %常规波束成型输出功率

s.P = diag(abs(s.P));                                      %提取特征值
s.P = s.P/max(s.P);                                        %输出功率归一化
s.P = 20*log10(abs(s.P));
end

%最小均方误差准则： 在非雷达应用中， 阵列协方差矩阵中通常都含有期望信
%号， 基于此种情况提出的准则。 使阵列输出与某期望响应的均方误差最小， 不需
%要知道期望信号的波达方向。 但是， 必须知道所需要的期望响应。
%最小均方差要求协方差中不含期望信号， 构造接收数据就是构造干扰加上噪声。 
function s = MMSE(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%数据处理%%%%%%%%%%%%%%%%%%%%%%%%%%
s.X = s.Xs +s.Xi+ s.Nos;                                   %构造接收数据
s.Rx = s.X*s.X'/s.L;                                       %接收数据相关矩阵
s.r_xd = s.X*s.ds';                                        %接收数据和期望信号互相关矩阵
Wopt = s.Rx'*s.r_xd;                                       %最佳权向量公式

for t = 1:length(s.Theta)
    a = exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*180/pi));
    s.P(t) = Wopt'*a;
end
s.P = abs(s.P);
s.P = s.P/max(max(s.P));
s.P = 20*log10(abs(s.P));
end

%最大信噪比准则： 使期望信号分量功率与噪声分量功率之比最大， 但是必须
%知道噪声的统计量和期望信号的波达方向。
function s=MSNR(s)
%%%%%%%%%%%数据处理%%%%%%%%%%%%
s.Rx=s.Xs*s.Xs'/s.L; %信号自相关矩阵
s.Rn=(s.Xi*s.Xi'+s.Nos*s.Nos')/s.L; %干扰+噪声自相关矩阵
[V,D]=eig(s.Rx,s.Rn); %广义特征值特征向量分解
[D,I]=sort(diag(D)); %排序Wopt=V(:,I(16)); %取最大特征值对应特征向量
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end


%线性约束最小方差准则： 对有用信号形式和来波方向完全已知， 在某种约束
%条件下使阵列输出的方差最小。
function s=LCMV(s)
%%%%%%%%%%%数据处理%%%%%%%%%%%%
s.X = s.Xi+s.Nos; %干扰+噪声
s.R = s.X*s.X'/s.L; %自相关协方差矩阵
%w=inv(R)*c/(c'*inv(R)*c)*f
Wopt = pinv(s.R)*s.As/(s.As'*pinv(s.R)*s.As); %最佳权向量
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end

%最大信干噪比准则： 使期望信号功率与干扰功率及噪声分量功率之和的比最大。
%必须知道干扰噪声统计数据和期望信号的波达方向。
function s=MSINR(s)
s.X = s.Xi+s.Nos; %干扰+噪声构造接收数据
s.R = s.X*s.X'/s.L; %构造干扰数据协方差矩阵
Wopt=pinv(s.R)*s.As; %MSINR最佳权向量算法
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end