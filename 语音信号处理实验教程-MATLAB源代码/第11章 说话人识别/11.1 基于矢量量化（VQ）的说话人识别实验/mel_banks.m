function H=mel_banks(p,N,fs,fl,fh)
%MEL_BANKS������mel�˲�����
% p:mel �˲�����   N��fft����  fs������Ƶ��
% H��Ϊ���ص��˲����飬ÿһ��Ϊһ���˲���
for m=1:p+2
    f(m)=(N/fs)*mel2freq(freq2mel(fl)+(m-1)*(freq2mel(fh)...
        -freq2mel(fl))/(p+1));
end
y=zeros(1,N);
H=zeros(p,N);
for m=2:p+1
    for k=1:N
        if f(m-1)<=k&&k<=f(m)
            y(k)=(k-f(m-1))/(f(m)-f(m-1));
        elseif f(m)<k&&k<=f(m+1)
            y(k)=(f(m+1)-k)/(f(m+1)-f(m));
        else
            y(k)=0;
        end  
    end
     H(m-1,:)=y;
end