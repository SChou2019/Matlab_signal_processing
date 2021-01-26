% x(n)���ο��ź�
% u(n)�������ź�
% d(n)�������ź�
% y(n)������ź�
% r(n)��x�˲�����ź�
% e(n)������ź�
% W(z)������Ӧ�˲���
% G(z)����ʵ�Ĵμ�ͨ��
% $\hat{G}(z)$�����ƵĴμ�ͨ��

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

% LMS��С�������
% S(z)�μ�ͨ�����ݺ���      % ys(k)�μ���Դ
% P(z)��ͨ�����ݺ���        % yp(k)������Դ
% C(z)������               % yw(k)������
% Sh(z)����������          % xs(k)�������ο��ź�

clear
T=1000; % �������ʱ��

% ���ǲ�֪��p(z)��S(z)���������Ǳ��뽨��dummy����·��
Pw=[0.01 0.25 0.5 1 0.5 0.25 0.01];
Sw=Pw*0.25;

x_iden=randn(1,T); % ����shape=(1,1000)�İ������źŹ���S(z)

% ����actuatorִ�У��ڴ�����λ�ò�����
y_iden=filter(Sw, 1, x_iden);

% Ȼ�󣬿�ʼʶ�����
Shx=zeros(1,16);       % ������Sh(z)��״̬
Shw=zeros(1,16);       % ������Sh(z)��Ȩ��
e_iden=zeros(1,  T);   % ʶ���������ݻ�����

%LMS �㷨
% [Shy,Shw]=lms(Shx,y_iden,x_iden,Shw,e_iden,T);
mu=0.1;                         % ѧϰ��
for k=1:T                      % ��ɢʱ�� k
    Shx=[x_iden(k) Shx(1:15)];  % ���´�������״̬
    Shy=sum(Shx.*Shw);            % ���㴫����Sh(z)�����
    e_iden(k)=y_iden(k)-Shy;    % �������     
    Shw=Shw+mu*e_iden(k)*Shx;   % ����Ȩ��
end

% �����
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
legend('S(z)ϵ��', 'Sh(z)ϵ��')


% ��second task������������������
X=randn(1,T);

% ����������λ�ý��յ�������
Yd=filter(Pw, 1, X);

% ����ϵͳ
Cx=zeros(1,16);       % C(z)��״̬
Cw=zeros(1,16);       % C(z)��Ȩ��
Sx=zeros(size(Sw));   % secondary��·��������״̬
e_cont=zeros(1,T);    % ���ƴ�������ݻ�����
Xhx=zeros(1,16);      % ���˺�x(k)��״̬

% FxLMS �㷨
% [Cy,Cw]=FxLMS(X,Cx,Cw,Sx,Sw,Shx,Shw,e_cont,Xhx,T,Yd);
mu=0.1;                            % ѧϰ��
for k=1:T                          % ��ɢʱ�� k
    Cx=[X(k) Cx(1:15)];            % ���¿�����״̬ 
    Cy=sum(Cx.*Cw);                % ������������
    Sx=[Cy Sx(1:length(Sx)-1)];    % ������secondary path
    e_cont(k)=Yd(k)-sum(Sx.*Sw);   % �����в�
    Shx=[X(k) Shx(1:15)];          % ����Sh(z)��״̬
    Xhx=[sum(Shx.*Shw) Xhx(1:15)]; % ������˺��x(k)
    Cw=Cw+mu*e_cont(k)*Xhx;        % ����controller��Ȩ��
end

% ���
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
legend('�����ź�', '�����ź�') 
