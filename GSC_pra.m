clc
clear all;
m = 8;%array阵元，形状类型？默认直线型
p = 4;%signal number信号数
N = 3000;%recursive number 迭代次数 或块拍数
A = zeros(m,p);%array apttern阵列流型
theta = [30 0 -45 60]*pi/180;% the signal from the direction of 30 degree is expected. DOA 30是期望信号方向
j = sqrt(-1);
w = [0.01 0.2 0.3 0.4]*pi; %frequency for each signal.各信号的数字频率，归一化表示

%% 
s = to_get_s(w,N,p);
s_rec = get_s_rec(s,m,p,theta);
S = s_rec;
%%% -------------------------%%自适应调节
y = s; % input data of GSC;
ad = exp(-j*pi*[0:m-1]'*sin(theta(1)));%steering vector in the direction expected. 期望信号方向导向矢量
c = 10;% conditio  波束成形条件
C = ad';
Wc = C'*inv(C*C')*C; %main path weight 主通道固定权重
Wa(1:m-1) = 0; %auxiliary path  辅助通道自适应权重
B = get_B(m,theta); %get Block Matrix 得到阻塞矩阵
u = 0.000001;
for k = 1:N
    yb = conj(B)*y(:,k);% m-1 * 1的列向量
    Zc(k) = Wc.'* y(:,k);% 转置的位置正确吗？？.不需要
    Za(k) = wa(:,k).'*yb;
    z(k) = Zc(k) - Za(k);
    wa(:,k+1) = wa(:,k) - u *z(k)*conj(yb);
end
%%%---------
%%%main path 主通道
wop = Wc;
drawpp(m,wop);
%%%auxiliary path 辅助通道
wop = B'*wa(:,N);
drawpp(m,wop);
%%array resopond 总的阵列响应
wop = Wc - B' * wa(:,N);
drawpp(m,wop);

function s = to_get_s(w,N,p)
s = zeros(p,N);
for i = 1:p
    s(i,1:N) = exp(1i*w(i).*(1:N));
end
end

function s_rec = get_s_rec(s,m,p,theta) %用于产生经过阵元后的信号数据
A = zeros(m,p);%  非必要代码
j = sqrt(-1);
%%% 阵元间距为半个波长
wi = pi*sin(theta);
A = exp(-j*wi'*[0:(m-1)]); %阵列流型
s_rec = A'*s;
s_rec = awgn(s_rec,10);
end

function y = drawpp(m,wop)
%figure,plot(e)?
thetas = [-90:90]; %角度间隔为1度
tm = thetas*pi/180;
am = exp(-j*pi*[0:m-1]'*sin(tm));
A = abs(wop.'*am);   %阵列响应
A = A/max(A);
figure,polar(tm,A)
A = 10*log10(A);   %对数图
hold on,title('归一化阵列响应极坐标图，八个阵元，信噪比20dB')
figure,plot(thetas,A);
hold on, title("八个阵元，信噪比20dB")
hold on, xlabel("入射角/度")
hold on,ylabel('归一化 A = 10*log10(A);')
grid on
axis ([-90 90 -35 0])
hold 
hold on,plot(-45,-35:0.1:0,'r');
hold on,plot(30,-35:0.1:0,'r');
hold on,plot(0,-35:0.1:0,'r');
hold on,plot(60,-35:0.1:0,'r');

end
