%��ʱƽ�����Ȳ��
function para=STAmdf(X)
    para=zeros(size(X));
    fn=size(X,2);              % ���֡��
    wlen=size(X,1);
    for i=1 : fn
        u=X(:,i);              % ȡ��һ֡
        for k=1:wlen
            para(:,k)=sum(abs(u(k:end)-u(1:end-k+1)));      %��ÿ��������ķ��Ȳ�
        end
    end
 end