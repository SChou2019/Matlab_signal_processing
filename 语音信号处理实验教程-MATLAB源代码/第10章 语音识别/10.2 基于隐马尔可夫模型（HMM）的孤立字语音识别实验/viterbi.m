function [prob,q]=viterbi(hmm,O)
%Viterbi�㷨
%���룺
%hmm--hmmģ��
%O--����۲����У�T*D��TΪ֡����DΪ����ά��
%�����
%prob--�������
%q--״̬����

init=hmm.init; %��ʼ����
trans=hmm.trans; %ת�Ƹ���
mix=hmm.mix;  %��˹���
N=hmm.N;  %HMM��״̬��
T=size(O,1); %����֡��

%����log(init)
ind1=find(init>0);
ind0=find(init<=0);
init(ind1)=log(init(ind1));
init(ind0)=-inf;

%����log(trans);
ind1 = find(trans>0);
ind0 = find(trans<=0);
trans(ind0) = -inf;
trans(ind1) = log(trans(ind1));

%��ʼ��
delta=zeros(T,N); %֡����״̬��
fai=zeros(T,N);
q=zeros(T,1);

%t=1;
for i=1:N
    delta(1,i)=init(i)+log(mixture(mix(i),O(1,:)));
end

%t=2:T
for t=2:T
    for j=1:N
        [delta(t,j),fai(t,j)]=max(delta(t-1,:)+trans(:,j)');
        delta(t,j)=delta(t,j)+log(mixture(mix(j),O(t,:)));
    end
end

%���ո��ʺ����ڵ�
[prob q(T)]=max(delta(T,:));

%�������״̬·��
for t=T-1:-1:1
    q(t)=fai(t+1,q(t+1));
end