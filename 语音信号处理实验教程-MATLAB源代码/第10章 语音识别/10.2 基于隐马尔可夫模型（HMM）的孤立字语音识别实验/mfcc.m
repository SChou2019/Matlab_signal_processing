function cc=mfcc(x)

%%-------׼������-------------
%��һ��mel�˲�����ϵ��(24����)
bank=melbankm2(24,256,8000,0,0.5,'m'); 
bank=full(bank);
bank=bank/max(bank(:));

%DCTϵ����12(�����mfcc����)��24
for k=1:12
 n=0:23;
 dctcoef(k,:)=cos(pi*k*(2*n+1)/(2*24));
end

%��һ���ĵ�����������
w=1+6*sin(pi*[1:12]./12);
w=w/max(w);

%--------��ȡ����-------------
%Ԥ�����˲���
xx=double(x);
xx=filter([1 -0.9375],1,xx);

%�����źŷ�֡
xx=enframe(x,256,80);

%����ÿ֡��MFCC����
for i=1:size(xx,1)
  y=xx(i,:);
  s=y'.*hamming(256);
  t=abs(fft(s));
  t=t.^2;
  c1=log(bank*t(1:129)); 
  c1=dctcoef*c1;
  c2=c1.*w';
  m(i,:)=c2';
end

%���ϵ��
dtm=zeros(size(m));
for i=3:size(m,1)-2
  dtm(i,:)=-2*m(i-2,:)-m(i-1,:)+m(i+1,:)+2*m(i+2,:);
end
dtm=dtm/3;

%�ϲ���֡����24��������������
cc=[m dtm];
cc=cc(3:size(m,1)-2,:);