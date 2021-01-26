clear all;close all;clc
%%%%%%%%%%%%��������%%%%%%%%%%%%
s.M = 16; %��Ԫ��
s.lamda = 1; %���ò���
s.d = 1/2*s.lamda; %ȡ1/2����Ϊ��Ԫ��
s.SNR = 10; %�����
s.INR = 10; %�����
s.Noise = 1; %��˹��������
s.SigDOA = [30 45]; %�ź�
s.IntDOA = [0 50 70]; %����
s.Theta = -90:1:90; %ɨ�跶Χ
s.L = 200; %������/������
s.fc = 100; %�ź�Ƶ��
s.fs=10000; %����Ƶ��
%%%%%%%%%%�����ź�%%%%%%%%%%%%
s.SigNum = length(s.SigDOA); %�ź���
s.IntNum = length(s.IntDOA); %������
s.ds = exp(1j*2*pi*(0:s.SigNum-1)'*s.fc/s.fs*(0:s.L-1)); %�ź�?
s.di = randn(s.IntNum,s.L)+1j*randn(s.IntNum,s.L); %����
s.As = exp(-1j*2*pi*s.d/s.lamda*(0:s.M-1)'*sin(s.SigDOA*pi/180)); %�źŷ���ʸ��
s.Ai = exp(-1j*2*pi*s.d/s.lamda*(0:s.M-1)'*sin(s.IntDOA*pi/180)); %���ŷ���ʸ��
s.Xs = sqrt(10^(s.SNR/10))*s.As*s.ds; %�����ź�
s.Xi = sqrt(10^(s.INR/10))*s.Ai*s.di; %�������
s.Nos = (randn(s.M,s.L)+1j*randn(s.M,s.L))/sqrt(2); %��������
s.Y = s.Xs+s.Xi+s.Nos; %�ϳ����н����ź�

%���沨���γ��ǶԲ������ŵ��źŽ��в����γɣ� ����������ݾ��ǹ����źż���������
%���沨���γ������ڵ����ź��Ҳ������ŵ���Դ�źŲ����γɡ� ����׼��˵
%�����沨���γɷ��Ĺ���ȱ����ǽǷֱ��ʵͣ� Ҫ��߽Ƿֱ��ʾ���Ҫ������Ԫ
%�����������Ԫ������
function s = CBF(s)
s.X = s.Xs + s.Nos;                                       %�����������
s.R = s.X*s.X'/s.L;                                       %Э�����
a = exp(-1i*pi*(0:s.M-1)*sin(s.Theta*pi/180));            %ɨ��
s.P = a'*s.R*a;                                           %���沨�������������

s.P = diag(abs(s.P));                                      %��ȡ����ֵ
s.P = s.P/max(s.P);                                        %������ʹ�һ��
s.P = 20*log10(abs(s.P));
end

%��С�������׼�� �ڷ��״�Ӧ���У� ����Э���������ͨ��������������
%�ţ� ���ڴ�����������׼�� ʹ���������ĳ������Ӧ�ľ��������С�� ����
%Ҫ֪�������źŵĲ��﷽�� ���ǣ� ����֪������Ҫ��������Ӧ��
%��С������Ҫ��Э�����в��������źţ� ����������ݾ��ǹ�����ż��������� 
function s = MMSE(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ݴ���%%%%%%%%%%%%%%%%%%%%%%%%%%
s.X = s.Xs +s.Xi+ s.Nos;                                   %�����������
s.Rx = s.X*s.X'/s.L;                                       %����������ؾ���
s.r_xd = s.X*s.ds';                                        %�������ݺ������źŻ���ؾ���
Wopt = s.Rx'*s.r_xd;                                       %���Ȩ������ʽ

for t = 1:length(s.Theta)
    a = exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*180/pi));
    s.P(t) = Wopt'*a;
end
s.P = abs(s.P);
s.P = s.P/max(max(s.P));
s.P = 20*log10(abs(s.P));
end

%��������׼�� ʹ�����źŷ���������������������֮����� ���Ǳ���
%֪��������ͳ�����������źŵĲ��﷽��
function s=MSNR(s)
%%%%%%%%%%%���ݴ���%%%%%%%%%%%%
s.Rx=s.Xs*s.Xs'/s.L; %�ź�����ؾ���
s.Rn=(s.Xi*s.Xi'+s.Nos*s.Nos')/s.L; %����+��������ؾ���
[V,D]=eig(s.Rx,s.Rn); %��������ֵ���������ֽ�
[D,I]=sort(diag(D)); %����Wopt=V(:,I(16)); %ȡ�������ֵ��Ӧ��������
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end


%����Լ����С����׼�� �������ź���ʽ������������ȫ��֪�� ��ĳ��Լ��
%������ʹ��������ķ�����С��
function s=LCMV(s)
%%%%%%%%%%%���ݴ���%%%%%%%%%%%%
s.X = s.Xi+s.Nos; %����+����
s.R = s.X*s.X'/s.L; %�����Э�������
%w=inv(R)*c/(c'*inv(R)*c)*f
Wopt = pinv(s.R)*s.As/(s.As'*pinv(s.R)*s.As); %���Ȩ����
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end

%����Ÿ����׼�� ʹ�����źŹ�������Ź��ʼ�������������֮�͵ı����
%����֪����������ͳ�����ݺ������źŵĲ��﷽��
function s=MSINR(s)
s.X = s.Xi+s.Nos; %����+���������������
s.R = s.X*s.X'/s.L; %�����������Э�������
Wopt=pinv(s.R)*s.As; %MSINR���Ȩ�����㷨
for t = 1:length(s.Theta)
a=exp(-1i*pi*(0:s.M-1)'*sin(s.Theta(t)*pi/180));
s.P(t) = Wopt'*a;
end
s.P=abs(s.P);
s.P=s.P/max(max(s.P));
s.P=20*log10(abs(s.P));
end