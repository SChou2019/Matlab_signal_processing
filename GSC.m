clc
clear all;
close all;

m=8 ;% array��Ԫ
p=4; %  signal number�ź���
N=3000;% recursive number �������� �������
A=zeros(m,p); % array pattern��������
theta=[30 0 -45 60]*pi/180;% the signal from the direction of 30 degree is expected. DOA 30Ϊ�����źŷ���
j=sqrt(-1);
w=[0.01 0.2 0.3 0.4]*pi; % frequency for each signal.�����źŵ�����Ƶ��

%
s=to_get_s(w,N,p);%����ԭʼ�ź�
s_rec=get_s_rec(s,m,p,theta);%���վ���
S=s_rec; %  output date matrix  .m*N ����Ԫ������ݾ���
%%%%����������������������������%% ����Ӧ����Ȩ
y=S;% input data of GSC;
ad=exp(-j*pi*[0:m-1]'*sin(theta(1))); %steering vector in the direction of expected. �����źŷ�����ʸ��
c=10;%  condition �����γ�����
%Ϊʲô��10���������ź�������Ԫ���й�ϵ��
C=ad';
Wc=C'*inv(C*C')*c; %main path weight ��ͨ���̶�Ȩ���Ի�������ź�
wa(1:m-1,1)=0;  % auxiliary path  ����ͨ������ӦȨ�������������͸����ź�
B=get_B(m,theta); % get Block Matrix �õ����������㷨��һ���о�
u=0.000001;
for k=1:N
    yb=conj(B)*y(:,k);  % m-1*1 ��������
    Zc(k)=Wc.'*y(:,k);
    Za(k)=wa(:,k).'*yb;
    Z(k)=Zc(k)-Za(k);
    wa(:,k+1)=wa(:,k)-u*Z(k)*conj(yb);
end
%%%%------------
%%%main path ��ͨ��
wop=Wc;
drawpp(m,wop);
%%%%auxiliary path ����ͨ��
wop=B'*wa(:,N);
drawpp(m,wop);
%%array response �ܵ�������Ӧ
wop=Wc-B'*wa(:,N);
drawpp(m,wop);


function y=drawpp(m,wop)
%figure,plot(e)
thetas=[-90:90];
tm=thetas*pi/180;
am=exp(-j*pi*[0:m-1]'*sin(tm));
A=abs(wop'*am);  %������Ӧ
A=A/max(A);
figure,polar(tm,A)
A=10*log10(A);  %����ͼ
hold on,title('��һ��������Ӧ��ֵ������ͼ������Ԫ�������20db')
figure,plot(thetas,A);
hold on,title('����Ԫ�������20db')
hold on,xlabel('�����/��')
hold on,ylabel('��һ�� A=10*log10(A);')
grid on 
axis([-90 90 -35 0]);
hold on,plot(-45,-35:0.1:0,'r');
hold on,plot(30,-35:0.1:0,'r');
hold on,plot(0,-35:0.1:0,'r');
hold on,plot(60,-35:0.1:0,'r');
end

function Bm=get_B(m,theta)  %���ڲ�����������%����������������������
u0=0.5*sin(theta(1)); % ������Ԫ���Ϊ�������
a0=exp(-j*2*pi*[0:m-1]'*u0);
u=u0+[1:m-1];
B=exp(-j*2*pi*[0:m-1]'*u);
Bm=conj(B');%% M-1*M �ľ���
end

function Rxx=get_Rxx(s_rec,N,p,m)
Rxx=zeros(m,m);
x=1:m;
for i=1:N
    for t=1:m
        x(t)=s_rec(t,i);
    end
    R=x'*x;
    Rxx=Rxx+R;
end
Rxx=Rxx/N;
end

function s_rec=get_s_rec(s,m,p,theta)  %���ڲ���������Ԫ����ź�����
A=zeros(m,p);
j=sqrt(-1);
%%% ��Ԫ���Ϊ�������
wi=pi*sin(theta);
A=exp(-j*wi'*[0:(m-1)]);  % ��������
s_rec=A'*s;
s_rec=awgn(s_rec,10);  % SNR=10 db
end

function s=to_get_s(w,N,p)
s=zeros(p,N);
for i=1:p
    s(i,1:N)=exp(j*w(i).*(1:N)); % ��ָ���ź�  �����ŵ�����Ϊ 1
end
end


