function mfc=my_mfcc(x,fs)
%MY_MFCC:��ȡSpeaker recognition�Ĳ���
%x:���������ź�  fs��������
%mfc��ʮ����mfccϵ����һ������ �Լ�һ�׺Ͷ��ײ�� ��36������ ÿһ��Ϊһ֡����
%clc,clear
%[x,fs]=audioread('wo6.wav');
N=256;p=24;
bank=mel_banks(p,N,fs,0,fs/2);
% ��һ��mel�˲�����ϵ��
bank=full(bank);
bank=bank/max(bank(:));

% ��һ��������������
w = 1 + 6 * sin(pi * [1:12] ./ 12);
w = w/max(w);

% Ԥ�����˲���
x=double(x);
x=filter([1 -0.9375],1,x);

% �����źŷ�֡
sf=check_ter(x,N,128,10);

x=div_frame(x,N,128);
%sf=sp_ter(x,4,16);
%m=zeros(size(x,1),13);
m=zeros(1,13);
% ����ÿ֡��MFCC����
for i=1:size(x,1)
    if sf(i)==1
        % j=j+1;
         y = x(i,:);
         y = y' .* hamming(N);
         energy=log(sum(y.^2)+eps);%����
         y = abs(fft(y));
         y = y.^2+eps;
        c1=dct(log(bank * y));
        c2 = c1(2:13).*w';%ȡ2~13��ϵ��
        %m(i,:)=[c2;energy]';
        m1=[c2;energy]';
        %m1=c2';
        m=[m;m1];
    end
    
end

%���ϵ��
dm = zeros(size(m));
dmm= zeros(size(m));
for i=2:size(m,1)-1
  dm(i,:) = (m(i,:) - m(i-1,:));
end
for i=3:size(m,1)-2
  dmm(i,:) = (dm(i,:) - dm(i-1,:));
end
%dm = dm / 3;

%�ϲ�mfcc������һ�ײ��mfcc����
mfc = [m dm dmm];
%ȥ����β��֡����Ϊ����֡��һ�ײ�ֲ���Ϊ0
mfc = mfc(3:size(m,1)-2,:);
%mfc=sum(mfc,1)/size(mfc,1);