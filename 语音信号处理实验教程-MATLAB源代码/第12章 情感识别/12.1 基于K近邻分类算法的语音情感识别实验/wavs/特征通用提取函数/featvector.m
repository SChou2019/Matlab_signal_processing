function feature=featvector(filename)
[y,fs]=wavread(filename); 
L=length(y);
ys=y;
for i=1:(length(y)-1)
    if (abs(y(i))<1e-3)  %  �޳���Сֵ�������ʱ����ʱʹ��  %
        ys(i)=ys(i+1);
        L=L-1;
    end
end 
y1=ys(1:L);
s=enframe(y,hamming(256),128); %  ��֡�Ӵ�  %
s1=enframe(y1,hamming(256),128); 
[nframe,framesize]=size(s);  
[nframe1,framesize1]=size(s1);
E=zeros(1,nframe1);  
Z=zeros(1,nframe);
F=zeros(1,nframe);
for i=1:nframe
    Z(i)=sum(abs(sign(s(i,framesize:2)-s(i,framesize-1:1))))/2;  %  ������  %
end
for i=1:nframe1
    E(i)=sum(s1(i,:).*s1(i,:)); %  ��ʱ����  %
end
%  ����Ƶ��  %
N=2048;R=4;
for i=1:nframe
    k = 1:R:N/2; K = length(k);  %  N��FFT�任������R�ǳ˵Ĵ�����f�ǲ���Ƶ��  %
    X = fft (s(i,:), N); 
    X=abs(X);  %  ��X������ֵ��ȡ������  %
    HPSx = X(k); 
    for r= R-1:-1:1
        HPSx = HPSx.*X (1:r:r*K);
    end
    [~,I]=max(HPSx);  %  ȡ���ֵ�㣬I�Ƕ�Ӧ�±�  %
    F(i)=I/N*fs; %  ����Ƶ��  %
end
%  ����֡��ֻ���  %
nf=1;
for i=1:(nframe-1)
    if(F(i)*F(i+1)~=0)
        dF(nf)=F(i)-F(i+1);
        nf=nf+1;
    end
end
%  0-250hz��ռ����  %
[s2,f1,t1]=specgram(y1,256,fs);
sn=20*log10(abs(s2)+eps);
sn1=sn+min(sn(:));
n=round(length(f1)*250/max(f1(:)));
Eratio=sum(sum(sn1(1:n,:)))/sum(sn1(:));

%{
%  Ƶ�׷���  %
out_gbvs=spGbvs(sn);  %  ������ͼ����ָ����H  %
H=out_gbvs.master_map_resized;
H(H<0.5)=0;
H(H>=0.5)=1;
snew=sn.*H;  %  �ָ��Ĳ���Ƶ��ͼ  %
row=size(snew,1);
col=size(snew,2);
ssnew=snew';
I=find(ssnew==0);
fa=ceil(I(1)/col);
fb=ceil(I(end)/col);
F1=f1(fb)-f1(fa);
%}

%  ���ƹ����  %
[fm,~] = formant_get(y,fs);
Fm1=fm(:,1);
Fm2=fm(:,2);
Fm3=fm(:,3);
%  MFCC  %
MFCCs=melcepst(y,fs,'0d'); %  MFCC����һ�ײ��ϵ��  %


%%    ������������    %%
%  ��ʱ����E  %
dim_max=140;
feature=zeros(dim_max,1);
x=0;t=0;
for i=1:(nframe1-1)
    t=abs(E(i)-E(i+1))/(nframe1-1);
    x=x+t;
end
E_shimmer=x/mean(E);
x1=0;x2=0;x3=0;x4=0;
for i=1:nframe1
    t1=i*mean(E);t2=i*E(i); t3=i*i;t4=i;
    x1=x1+t1;x2=x2+t2;x3=x3+t3;x4=x4+t4;
end
x4=x4*x4/nframe1;
s1=x2-x1;s2=x3-x4;
E_Reg_coff=s1/s2;
x=0;
for i=1:nframe1
    t=E(i)-(mean(E)-s1/s2*x4/nframe1)-s1/s2*i;
    x=x+t^2/nframe1;
end
E_Sqr_Err=x;
feature(1:7,1)=[max(E);min(E);mean(E);var(E);E_shimmer;E_Reg_coff;E_Sqr_Err];%  ��ʱ�����������  %

%  ������  %
feature(8,1)=Eratio;

%  ����Ƶ��F  %
x=0;
for i=1:(nframe-1)
    t=abs(F(i)-F(i+1));
    x=x+t;
end
F_Jitter1=100*x/(mean(F)*(nframe-1));
x=0;
for i=2:(nframe-1)
    t=abs(2*F(i)-F(i+1)-F(i-1));
    x=x+t;
end
F_Jitter2=100*x/(mean(F)*(nframe-2));

%% ʹF����Сֵ����Ч��ȥ����ֵ��
k=1;
for i=2:numel(F)
    if(F(i)==F(1))
        continue;
    end
    FF(k)= F(i);
    k=k+1;
 
end

feature(9:14,1)=[max(F);min(FF);mean(F);var(F);F_Jitter1;F_Jitter2];%  ����Ƶ���������  %

%  ����֡��ֻ���  %
feature(15:18,1)=[max(dF);min(dF);mean(dF);var(dF)];%  ����֡��ֻ���  %

%  �����  %
x1=0;x2=0;x3=0;
for i=1:(numel(Fm1)-1)
    t1=abs(Fm1(i)-Fm1(i+1));
    t2=abs(Fm2(i)-Fm2(i+1));
    t3=abs(Fm3(i)-Fm3(i+1));
    x1=x1+t1;x2=x2+t2;x3=x3+t3;
end
Fm1_Jitter1=100*x1/(mean(Fm1)*(numel(Fm1)-1));%  ǰ����������һ�׶���  %
Fm2_Jitter1=100*x2/(mean(Fm2)*(numel(Fm1)-1));
Fm3_Jitter1=100*x3/(mean(Fm2)*(numel(Fm1)-1));
Fm2R=Fm2./(Fm2-Fm1);
nFm=[max(Fm1);min(Fm1);mean(Fm1);var(Fm1);Fm1_Jitter1;max(Fm2);min(Fm2);mean(Fm2);var(Fm2);Fm2_Jitter1;max(Fm3);min(Fm3);mean(Fm3);var(Fm3);Fm3_Jitter1;max(Fm2R);min(Fm2R);mean(Fm2R)];%  ������������  %
feature(19:(size(nFm)+18),1)=nFm;%  20-37  %
% size(feature)
%  MFCCs & dMFCCs  %
for i=1:size(MFCCs,2)
    feature(37+4*(i-1):37+4*i-1,1)=[max(MFCCs(:,i));min(MFCCs(:,i));mean(MFCCs(:,i));var(MFCCs(:,i))];%  mel����ϵ������һ�ײ���������  %
end




