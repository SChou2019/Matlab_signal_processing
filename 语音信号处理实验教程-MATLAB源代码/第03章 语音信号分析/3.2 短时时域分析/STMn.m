%��ʱƽ�����ȼ��㺯��
function para=STMn(x,win,inc)
    X=enframe(x,win,inc)';     % ��֡
    fn=size(X,2);              % ���֡��
    for i=1 : fn
        u=X(:,i);              % ȡ��һ֡
        para(i)=sum(abs(u))/200;         % ��һ֡�ۼ����
    end
 end