function H=melbankm(M,N,fs,fl,fh)
%Mel�˲�����
% �������:	M �˲���������
%		   N  FFT����
%		   fs  ����Ƶ��
%		   fl -fh Ϊ���Թ����׵�����Ƶ��(Ĭ��Ϊ0-0.5*fs)
%		   
%�������:	H�����˲����飬ÿһ��Ϊһ���˲���
if nargin<4
    fl=0*fs;
    fh=0.5*fs;
end

%����ÿ���˲���������Ƶ��
f=zeros(1,M+1);
for m=1:M+2
    f(m)=floor((N/fs)*mel2freq(freq2mel(fl)...
        +(m-1)*(freq2mel(fh)-freq2mel(fl))/(M+1)));
end
%���˲�����H
c=floor(N/2)+1;
y=zeros(1,c);
H=zeros(M,c);
for m=2:M+1
    for k=1:c                   %����fh���Ϊfs/2����ô�����cλ���ܴ洢H
        if f(m-1)<=k&&k<=f(m)
            y(k)=(k-f(m-1))/(f(m)-f(m-1));
        elseif f(m)<=k&&k<=f(m+1)
            y(k)=(f(m+1)-k)/(f(m+1)-f(m));
        else
            y(k)=0;
        end
    end
    H(m,:)=y;
end
    