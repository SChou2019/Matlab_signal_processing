%��ʱ�������㺯��
function para=STEn(x,win,inc)
    X=enframe(x,win,inc)';     % ��֡
    fn=size(X,2);              % ���֡��
    for i=1 : fn
        u=X(:,i);              % ȡ��һ֡
        u2=u.*u;               % �������
        para(i)=sum(u2);         % ��һ֡�ۼ����
    end
 end