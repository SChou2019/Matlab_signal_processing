%��ʱ�����ʼ��㺯��
function para=STZcr(x,win,inc)
    X=enframe(x,win,inc)';        % ��֡
    fn=size(X,2);                       % ���֡��
    wlen=length(win);               % ���֡��
    para=zeros(1,fn);                 % ��ʼ��
    for i=1:fn
        z=X(:,i);                           % ȡ��һ֡����
        for j=1: (wlen- 1) ;            % ��һ֡��Ѱ�ҹ����
             if z(j)* z(j+1)< 0         % �ж��Ƿ�Ϊ�����
                 para(i)=para(i)+1;   % �ǹ���㣬��¼1��
             end
        end
    end
 end