%��ʱ����غ���
function para=STAc(X)
    para=zeros(size(X));
    fn=size(X,2);                            % ���֡��
    wlen=size(X,1);                         %��֡��
    for i=1 : fn
        u=X(:,i);                               % ȡ��һ֡
        R=xcorr(u);                         %��ʱ����ؼ���
        para(:,i)=R(wlen,end);          %ֻȡkΪ��ֵ������غ���
    end
 end