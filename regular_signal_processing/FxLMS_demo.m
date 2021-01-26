% x(n)：参考信号
% u(n)：控制信号
% d(n)：期望信号
% y(n)：输出信号
% r(n)：x滤波后的信号
% e(n)：误差信号
% W(z)：自适应滤波器
% G(z)：真实的次级通道
% $\hat{G}(z)$：估计的次级通道

%              +-----------+                       +   
% x(k) ---+--->|   P(z)    |--yp(k)----------------> sum --+---> e(k)
%         |    +-----------+                          ^-   |
%         |                                           |    |
%         |        \                                ys(k)  |     
%         |    +-----------+          +-----------+   |    |
%         +--->|   C(z)    |--yw(k)-->|   S(z)    |---+    |
%         |    +-----------+          +-----------+        |
%         |            \                                   |
%         |             \----------------\                 |
%         |                               \                |
%         |    +-----------+          +-----------+        |
%         +--->|   Sh(z)   |--xs(k)-->|    LMS    |<-------+
%              +-----------+          +-----------+        

% LMS最小均方误差
% S(z)次级通道传递函数      % ys(k)次级声源
% P(z)主通道传递函数        % yp(k)初级声源
% C(z)控制器               % yw(k)控制器
% Sh(z)传感器函数          % xs(k)传感器参考信号

clear
T=1000; % 仿真持续时间

% 我们不知道p(z)和S(z)，所以我们必须建立dummy虚拟路径
Pw=[0.01 0.25 0.5 1 0.5 0.25 0.01];
Sw=Pw*0.25;

x_iden=randn(1,T); % 产生shape=(1,1000)的白噪声信号估计S(z)

% 送至actuator执行，在传感器位置测量，
y_iden=filter(Sw, 1, x_iden);

% 然后，开始识别过程
Shx=zeros(1,16);       % 传感器Sh(z)的状态
Shw=zeros(1,16);       % 传感器Sh(z)的权重
e_iden=zeros(1,  T);   % 识别错误的数据缓冲区

%LMS 算法
% [Shy,Shw]=lms(Shx,y_iden,x_iden,Shw,e_iden,T);
mu=0.1;                         % 学习率
for k=1:T                      % 离散时间 k
    Shx=[x_iden(k) Shx(1:15)];  % 更新传感器的状态
    Shy=sum(Shx.*Shw);            % 计算传感器Sh(z)的输出
    e_iden(k)=y_iden(k)-Shy;    % 计算误差     
    Shw=Shw+mu*e_iden(k)*Shx;   % 调整权重
end

% 检查结果
subplot(2,1,1)
plot((1:T), e_iden)
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Identification error');
subplot(2,1,2)
stem(Sw) 
hold on 
stem(Shw, 'r*')
ylabel('Amplitude');
xlabel('Numbering of filter tap');
legend('S(z)系数', 'Sh(z)系数')


% 第second task二个任务是主动控制
X=randn(1,T);

% 测量传感器位置接收的噪声，
Yd=filter(Pw, 1, X);

% 启动系统
Cx=zeros(1,16);       % C(z)的状态
Cw=zeros(1,16);       % C(z)的权重
Sx=zeros(size(Sw));   % secondary次路径的虚拟状态
e_cont=zeros(1,T);    % 控制错误的数据缓冲区
Xhx=zeros(1,16);      % 过滤后x(k)的状态

% FxLMS 算法
% [Cy,Cw]=FxLMS(X,Cx,Cw,Sx,Sw,Shx,Shw,e_cont,Xhx,T,Yd);
mu=0.1;                            % 学习率
for k=1:T                          % 离散时间 k
    Cx=[X(k) Cx(1:15)];            % 更新控制器状态 
    Cy=sum(Cx.*Cw);                % 计算控制器输出
    Sx=[Cy Sx(1:length(Sx)-1)];    % 传播到secondary path
    e_cont(k)=Yd(k)-sum(Sx.*Sw);   % 测量残差
    Shx=[X(k) Shx(1:15)];          % 更新Sh(z)的状态
    Xhx=[sum(Shx.*Shw) Xhx(1:15)]; % 计算过滤后的x(k)
    Cw=Cw+mu*e_cont(k)*Xhx;        % 调整controller的权重
end

% 结果
figure
subplot(2,1,1)
plot((1:T), e_cont)
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')
subplot(2,1,2)
plot((1:T), Yd) 
hold on 
plot((1:T), Yd-e_cont, 'r:')
ylabel('Amplitude');
xlabel('Discrete time k');
legend('噪声信号', '控制信号') 
