function hmm = inithmm(obs, N, M)

K=length(obs);  %���������������˴�Ϊ10��
% N=length(M)        %״̬�����˴�Ϊ4��
hmm.N=N;
hmm.M=M;           %ÿ��״̬��pdf��

%-------��ʼ����ʼ״̬����pai
hmm.init=zeros(N,1);
hmm.init(1)=1;

%-------��ʼ��״̬ת�Ƹ��ʾ���A
hmm.trans=zeros(N,N);
for i=1:N-1
  hmm.trans(i,i)=0.5;
  hmm.trans(i,i+1)=0.5;
end
hmm.trans(N,N)=1;

% -------��ʼ������۲�ֵ����B�����������̬�ֲ���
for k=1:K
 T=size(obs(k).fea,1);  
 obs(k).segment=floor([1:T/N:T T+1]);
end

for i=1:N
   vector=[];
   for k=1:K
       seg1=obs(k).segment(i);
       seg2=obs(k).segment(i+1)-1;
       vector=[vector;obs(k).fea(seg1:seg2,:)];
   end
   mix(i)=getmix(vector,M(i)); %����getmix���������ؽṹ��mix(i)
end

hmm.mix=mix;

%------------getmix������ʵ�ֹ���------
function mix=getmix(vector,M)

%����ÿ���ɷֵĳ�ʼ��ֵmiu
[mean esp nn]=kmeans1(vector,M);

% ����ÿ���ɷֵ�Э����sigma
for j=1:M
    ind=find(j==nn);
    tmp=vector(ind,:);
    var(j,:)=std(tmp);
end

% ����ÿ���ɷֵ�Ȩ��w
weight=zeros(M,1);
for j=1:M
  for k=1:size(nn,1)  
    if nn(k)==j
        weight(j)=weight(j)+1;
    end
  end
end
weight=weight/sum(weight);

% ������
mix.M=M;
mix.mean=mean;    
mix.var=var.^2;	
mix.weight=weight;	

