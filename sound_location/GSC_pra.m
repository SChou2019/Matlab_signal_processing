clc
clear all;
m = 8;%array��Ԫ����״���ͣ�Ĭ��ֱ����
p = 4;%signal number�ź���
N = 3000;%recursive number �������� �������
A = zeros(m,p);%array apttern��������
theta = [30 0 -45 60]*pi/180;% the signal from the direction of 30 degree is expected. DOA 30�������źŷ���
j = sqrt(-1);
w = [0.01 0.2 0.3 0.4]*pi; %frequency for each signal.���źŵ�����Ƶ�ʣ���һ����ʾ

%% 
s = to_get_s(w,N,p);
s_rec = get_s_rec(s,m,p,theta);
S = s_rec;
%%% -------------------------%%����Ӧ����
y = s; % input data of GSC;
ad = exp(-j*pi*[0:m-1]'*sin(theta(1)));%steering vector in the direction expected. �����źŷ�����ʸ��
c = 10;% conditio  ������������
C = ad';
Wc = C'*inv(C*C')*C; %main path weight ��ͨ���̶�Ȩ��
Wa(1:m-1) = 0; %auxiliary path  ����ͨ������ӦȨ��
B = get_B(m,theta); %get Block Matrix �õ���������
u = 0.000001;
for k = 1:N
    yb = conj(B)*y(:,k);% m-1 * 1��������
    Zc(k) = Wc.'* y(:,k);% ת�õ�λ����ȷ�𣿣�.����Ҫ
    Za(k) = wa(:,k).'*yb;
    z(k) = Zc(k) - Za(k);
    wa(:,k+1) = wa(:,k) - u *z(k)*conj(yb);
end
%%%---------
%%%main path ��ͨ��
wop = Wc;
drawpp(m,wop);
%%%auxiliary path ����ͨ��
wop = B'*wa(:,N);
drawpp(m,wop);
%%array resopond �ܵ�������Ӧ
wop = Wc - B' * wa(:,N);
drawpp(m,wop);

function s = to_get_s(w,N,p)
s = zeros(p,N);
for i = 1:p
    s(i,1:N) = exp(1i*w(i).*(1:N));
end
end

function s_rec = get_s_rec(s,m,p,theta) %���ڲ���������Ԫ����ź�����
A = zeros(m,p);%  �Ǳ�Ҫ����
j = sqrt(-1);
%%% ��Ԫ���Ϊ�������
wi = pi*sin(theta);
A = exp(-j*wi'*[0:(m-1)]); %��������
s_rec = A'*s;
s_rec = awgn(s_rec,10);
end

function y = drawpp(m,wop)
%figure,plot(e)?
thetas = [-90:90]; %�Ƕȼ��Ϊ1��
tm = thetas*pi/180;
am = exp(-j*pi*[0:m-1]'*sin(tm));
A = abs(wop.'*am);   %������Ӧ
A = A/max(A);
figure,polar(tm,A)
A = 10*log10(A);   %����ͼ
hold on,title('��һ��������Ӧ������ͼ���˸���Ԫ�������20dB')
figure,plot(thetas,A);
hold on, title("�˸���Ԫ�������20dB")
hold on, xlabel("�����/��")
hold on,ylabel('��һ�� A = 10*log10(A);')
grid on
axis ([-90 90 -35 0])
hold 
hold on,plot(-45,-35:0.1:0,'r');
hold on,plot(30,-35:0.1:0,'r');
hold on,plot(0,-35:0.1:0,'r');
hold on,plot(60,-35:0.1:0,'r');

end
