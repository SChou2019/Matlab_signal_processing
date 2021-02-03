function FMatrix=mfccf(num,s,Fs)
%���㲢�����ź�s��mfcc��������һ�׺Ͷ��ײ�ֲ���
%����˵���� num ��mfccϵ������FsΪ����Ƶ��

N=512;              % FFT ��
Tf=0.02;            %���ڵ�ʱ��
n=Fs*Tf;            %ÿ�����ڵĳ���
M=24;               %MΪ�˲�������
l=length(s);        
Ts=0.01;            %֡��ʱ��
FrameStep=Fs*Ts;    %֡��
a=1;
b=[1, -0.97];       %Ԥ���ش����һ��FIR��ͨ�˲���ϵ��


noFrames=fix((l-n)/FrameStep)+1;  %֡��
FMatrix=zeros(noFrames-2, num); %����ϵ����ʼ����һ�ײ�ֵ���β2֡Ϊ0,
lifter=1:num;                   %���׼�Ȩ��ʼ����
lifter=1+floor(num/2)*(sin(lifter*pi/num));%���׼�Ȩ����

if mean(abs(s)) > 0.01
    s=s/max(s);                     %��ʼ��
end

%���� MFCC ϵ��
for i=1:noFrames
    frame=s((i-1)*FrameStep+1:(i-1)*FrameStep+n);  %frame���ڴ洢ÿһ֡����
    framef=filter(b,a,frame);   %Ԥ�����˲���
    F=framef.*hamming(n);       %�Ӻ�����
    FFTo=fft(F,n);         %����FFT
    melf=melbankm(M,N,Fs);      %����mel�˲�����ϵ��  
    halfn=1+floor(N/2);    
    spectr=log(melf*(abs(FFTo(1:halfn)).^2)+1e-22);%�����˲��������˲�����1e-22��ֹ���ж��������������
%%%%%%%%%%����DCT�任%%%%%%%
c=zeros(1,num);
    for p=1:num
        for m=1:M
           c(p)=c(p)+spectr(m)*cos(p*(m-0.5)*pi/M);
        end
    end
    ncoeffs=c.*lifter;          %��MFCC�������е��׼�Ȩ
    FMatrix(i, :)=ncoeffs;     
end

%����deltacoeff��������MFCC���ϵ��
d=deltacoeff(FMatrix);         %����һ�ײ��ϵ��
d1=deltacoeff(d);              %������ײ��ϵ��
FMatrix=[FMatrix,d,d1];        %������������Ϊ��������
