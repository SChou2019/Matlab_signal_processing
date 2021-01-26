%����ź�DOA�����㷨��������ź��ӿռ䣨ISM���㷨 
% Developed by xiaofei zhang (�Ͼ����պ����ѧ ���ӹ���ϵ ��С�ɣ�
% EMAIL:zhangxiaofei@nuaa.edu.cn


clc 
clear all 
close all 
M=12;                                      %��Ԫ�� 
N=200;                                    %������ 
ts=0.01;                                  %ʱ�������� 
f0=100;                                   %�����ź�����Ƶ��  
f1=80;                                    %�����ź����Ƶ�� 
f2=120;                                   %�����ź����Ƶ�� 
c=1500;                                   %���� 
lambda=c/f0;                              %���� 
d=lambda/2;                               %��Ԫ��� 
SNR=15;                                   %����� 
b=pi/180;   
theat1=30*b;                               %�����źŲ�����1 
theat2=0*b;                                %�����źŲ�����2 
n=ts:ts:N*ts;                             
theat=[theat1 theat2]'; 
 
%%%%%%%%%%%%%%%% produce signal %%%%%%%%%%%%%%%% 
 
s1=chirp(n,80,1,120);                     %�������Ե�Ƶ�ź�1�� 
sa=fft(s1,2048);                          %����FFT�任 
%figure, %specgram(s1,256,1E3,256,250);   %Ƶ��ͼ
s2=chirp(n+0.100,80,1,120);                %�������Ե�Ƶ�ź�2 
sb=fft(s2,2048);                           %����FFT�任 
  
% %%%%%%%%%%%%%%%%%%%%% ISM�㷨 %%%%%%%%%%%%%%%%%% 
P=1:2; 
a=zeros(M,2); 
sump=zeros(1,181); 
for i=1:N
    f=80+(i-1)*1.0; 
    s=[sa(i) sb(i)]'; 
    for m=1:M 
        a(m,P)=exp(-j*2*pi*f*d/c*sin(theat(P))*(m-1))'; 
    end 
    R=a*(s*s')*a'; 
    [em,zm]=eig(R); 
    [zm1,pos1]=max(zm); 
    for l=1:2 
        [zm2,pos2]=max(zm1); 
        zm1(:,pos2)=[]; 
        em(:,pos2)=[]; 
    end 
    k=1; 
    for ii=-90:1:90 
        arfa=sin(ii*b)*d/c; 
        for iii=1:M 
            tao(1,iii)=(iii-1)*arfa; 
        end 
        A=[exp(-j*2*pi*f*tao)]'; 
        p(k)=A'*em*em'*A; 
        k=k+1; 
    end 
    sump=sump+abs(p); 
end 
pmusic=1/33*sump; 
pm=1./pmusic; 
thetaesti=-90:1:90; 
plot(thetaesti,20*log(abs(pm))); 
xlabel('�����/��'); 
ylabel('�ռ���/dB'); 
grid on
