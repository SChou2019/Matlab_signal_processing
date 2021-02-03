function hmm=baum_welch(hmm,obs)

mix=hmm.mix; %��˹���ģ��
N=hmm.N; %HMM��״̬��
K=length(obs); %ѵ������������
SIZE=size(obs(1).fea,2); %����ʸ���ĸ���

for loop = 1:40
   % ----����ǰ��, ������ʾ���
   for k=1:K
     param(k)=getparam(hmm,obs(k).fea);
   end

   %----�ع�ת�Ƹ��ʾ���A
   for i=1:N-1
     demon=0;
     for k=1:K  
        tmp=param(k).ksai(:,i,:);
        demon=demon+sum(tmp(:)); %��ʱ��t��j���
     end
     for j=i:i+1  
        nom=0;
        for k=1:K  
            tmp=param(k).ksai(:,i,j);
            nom=nom+sum(tmp(:));  %��ʱ��t���
        end
        hmm.trans(i,j)=nom/demon;
     end
   end

   %----�ع�����۲�ֵ����B
   for j=1:N %״̬ѭ��
     for l=1:hmm.M(j) %��ϸ�˹����Ŀ
        %�������ϳɷֵľ�ֵ��Э�������
        nommean=zeros(1,SIZE);
        nomvar=zeros(1,SIZE);
        denom=0;
        for k=1:K  %ѵ����Ŀ��ѭ��
           T=size(obs(k).fea,1);  %֡��
           for t=1:T   %֡����ʱ�䣩�ı���
             x=obs(k).fea(t,:);
             nommean=nommean+param(k).gama(t,j,l)*x;
             nomvar=nomvar+param(k).gama(t,j,l)*(x-mix(j).mean(l,:)).^2;
             denom=denom+param(k).gama(t,j,l);
           end
        end
        hmm.mix(j).mean(l,:)=nommean/denom;
        hmm.mix(j).var(l,:)=nomvar/denom;
   
        %�������ϳɷֵ�Ȩ��
        nom=0;
        denom=0;
        for k=1:K
          tmp=param(k).gama(:,j,l);
          nom=nom+sum(tmp(:));
          tmp=param(k).gama(:,j,:);
          denom=denom+sum(tmp(:));
        end
        hmm.mix(j).weight(l)=nom/denom;
     end
   end   
      
end
  




