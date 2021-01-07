% sound loaction
clc
clear all;
close all
%初始化设置
c = 334;%实际与温度相关
fs = 1000;%采样频率
T = 0.1; 
t = 0:1/fs:T;  % 时间 [0, 0.1]
L = length(t); % 时间长度:101
f = 500; %感兴趣的频率
w = 2*pi*f; %转换成角频率
k = w/c;%波长k
% 振元位置
M = 17;%阵元个数
% Nmid = 12;      % 参考点
% d = 3;         % 阵元间距
% m = (0:1:M-1) 
yi = zeros(M,1);%阵元均位于y平面上
zi = [ 0; 3; 6; 9;12;15;18;21;24;12;12;12;12;12;12;12;12];
xi = [12;12;12;12;12;12;12;12;12; 0; 3; 6; 9;15;18;21;24];
%xi = xi.'      % 列向量 m*d 阵元数*阵元间距
figure(1)
plot(xi,zi,'r*');
title('十字形麦克风阵列')

%声源位置，实际应用中不知道，扫描角度和间距确定

x1 = 12;
y1 = 10;
z1 = 12; %声源位置 （12,10,12） x,z为水平面，可以由俯、仰角确定
 
x2 = 12;  % array center
y2 = 0;
z2 = 12;
Ric1 = sqrt((x1-xi).^2+(y1-yi).^2+(z1-zi).^2);%声源到各阵元的距离
Ric2 = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2); %sound to array center:10
Rn1 = Ric1 - Ric2; %声源至各阵元与参考阵元的声程差矢量


% Ric1 = sqrt((x1-xi).^2+(y1-yi).^2+(z1-zi).^2); % 声源到各阵元的距离
% Ric2 = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2); %sound to array center:10
% Rn1 = Ric1 - Ric2; %声源至各阵元与参考阵元的声程差矢量

s1 = cos(2*w*t);%参卡很远接收到的矢量
Am = 10^(-1);%振幅
n1 = Am * (randn(M, L) + 1i*randn(M, L)); % 各阵元高斯白噪声
p1 = zeros(M,L);
% 各阵元的延时求和，同时球的各阵元接收的声压矩阵信号，以及协方差矩阵
for k1 = 1:M
    p1(k1,:) =Ric2/Ric1(k1)*s1.*exp(-1i*w*Rn1(k1)/c);
    %p1(k1,:) = Ric2/Ric1(k1) * s1.*exp(-j*w*Rn1(k1)/c);
end


p = p1 + n1;%加上噪声
R = p*p'/L;%接收数据的自协方差


%扫面范围
%步长设置0.1，扫描范围20*20平面，双重循环得到M*1矢量矩阵，最后得到交叉谱矩阵（cross spectrum matrix）
%由DSP理论，得到声功率
step_x = 0.1;
step_z = 0.1;
y = y1;
x = (9:step_x:15);
z = (9:step_z:15);
for k1 = 1:length(z)
    for k2 = 1:length(x)
        Ri = sqrt((x(k2)-xi).^2+(y-yi).^2+(z(k1)-zi).^2);  % 该扫描点到各阵元的聚焦距离矢量
        Ri2 = sqrt((x(k2)-x2).^2+(y-y2).^2+(z(k1)-z2).^2);  % 10.8628
        Rn = Ri - Ri2;%量程差
        b = exp(-1i*w*Rn/c);%声聚焦方向矢量
        Pcbf(k1,k2) = abs(b'*R*b);% ???CSM,最关键,(1,18)*(18,18)*(18,1)
    end
end



%归一化处理
for k1 = 1:length(z)
    pp(k1) = max(Pcbf(k1,:));%Pcbf每行最大值
end
Pcbf = Pcbf/max(pp);%全局归一化
figure(2)
surf(x,z,Pcbf)
xlabel('x(m)');zlabel('z(m)')
title('三维单声源图')
colorbar
 
figure(3)
pcolor(x,z,Pcbf);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
title('单声源图')
colorbar


