[x,fs]=audioread('C:\Users\lenovo\Desktop\�����ļ�\Ready.wav');
bank=melbankm(24,256,fs,0,0.4,'t');%Mel�˲����Ľ���Ϊ24��fft�任�ĳ���Ϊ256������Ƶ��Ϊ16000Hz
%��һ��mel�˲�����ϵ��
bank=full(bank);%ϡ�����ת��Ϊ��ɾ���
bank=bank/max(bank(:));
for k=1:12			%��һ��mel�˲�����ϵ��
    n=0:23;
    dctcoef(k,:)=cos((2*n+1)*k*pi/(2*24));
end
w=1+6*sin(pi*[1:12]./12);%��һ��������������
w=w/max(w);%Ԥ�����˲���
xx=double(x);
xx=filter([1 -0.9375],1,xx);%�����źŷ�֡,Ԥ���ز���
xx=enframe(xx,256,80);%��x 256���Ϊһ֡
%����ÿ֡��MFCC����
for i=1:size(xx,1)
    y=xx(i,:);
    s=y'.*hamming(256);
	t=abs(fft(s));%fft���ٸ���Ҷ�任
    t=t.^2;
	c1=dctcoef*log(bank*t(1:129));
    c2=c1.*w';
    m(i,:)=c2';
end
%��ȡһ�ײ��ϵ��
dtm=zeros(size(m));
for i=3:size(m,1)-2
    dtm(i,:)=-2*m(i-2,:)-m(i-1,:)+m(i+1,:)+2*m(i+2,:);
end
dtm=dtm/3;
%��ȡ���ײ��ϵ��
dtmm=zeros(size(dtm));
for i=3:size(dtm,1)-2
    dtmm(i,:)=-2*dtm(i-2,:)-dtm(i-1,:)+dtm(i+1,:)+2*dtm(i+2,:);
end
dtmm=dtmm/3;
%�ϲ�mfcc������һ�ײ��mfcc����
ccc=[m dtm dtmm];
%ȥ����β��֡����Ϊ����֡��һ�ײ�ֲ���Ϊ0
ccc=ccc(3:size(m,1)-2,:);
ccc;
subplot(2,1,1) 
ccc_1=ccc(:,1);
plot(ccc_1);title('MFCC');ylabel('��ֵ');
[h,w]=size(ccc);
A=size(ccc);
subplot(212)    
plot([1,w],A);
xlabel('ά��');
ylabel('��ֵ');
title('ά�����ֵ�Ĺ�ϵ') 
