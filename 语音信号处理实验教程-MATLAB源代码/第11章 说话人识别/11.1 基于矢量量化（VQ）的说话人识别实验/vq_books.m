clc,clear
%���˵����ʶ���ѵ����ƥ��
k=8;
N=4;%N��˵����
for i=1:N
    s=['SX',num2str(i),'.wav'];
    [x,fs]=audioread(s);
    x=x/max(x);%��һ��
    mel=my_mfcc(x,fs)';%ÿ��Ϊһ������
    v=lbg(mel,k);
    u{i}=[v(1:k).mea];
end
save data.mat u