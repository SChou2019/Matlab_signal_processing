%��ʱ�����ʼ��㺯��
function zcr=STZcr(x,win,inc)
    X=enframe(x,win,inc)';        % ��֡
    fn=size(X,2);                       % ���֡��
    if length(win)==1
        wlen=win;               % ���֡��
    else
        wlen=length(win);
    end
    zcr=zeros(1,fn);                 % ��ʼ��
    delta=0.01;                                % ����һ����С����ֵ
    for i=1:fn
        z=X(:,i);                           % ȡ��һ֡����
    for k=1 : wlen                         % ���Ľط�����
        if z(k)>=delta
            ym(k)=z(k)-delta;
        elseif z(k)<-delta
            ym(k)=z(k)+delta;
        else
            ym(k)=0;
        end
    end
    zcr(i)=sum(ym(1:end-1).*ym(2:end)<0);  % ȡ�ô�����һ֡����Ѱ�ҹ�����
    end
