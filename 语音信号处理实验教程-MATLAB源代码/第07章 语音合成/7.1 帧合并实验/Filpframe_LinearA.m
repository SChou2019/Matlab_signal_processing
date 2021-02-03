%���ڱ����ص���ӷ����źŻ�ԭ����
function frameout=Filpframe_LinearA(x,win,inc)

[nf,len]=size(x);
nx=(nf-1) *inc+len;                 %ԭ�źų���
frameout=zeros(nx,1);
nwin=length(win);                   % ȡ����
overlap=nwin-inc;                         % �ص�����
tempr1=(0:overlap-1)'/overlap;            % б���Ǵ�����w1
tempr2=(overlap-1:-1:0)'/overlap;         % б���Ǵ�����w2
if (nwin ~= 1)                           % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
    winx=repmat(win',nf,1);
    x=x./winx;                          % ��ȥ�Ӵ���Ӱ��
    x(find(isinf(x)))=0;                %ȥ����0�õ���Inf
end

for i=1:nf
    xn=x(i,:)';
    if i==1                           % ��Ϊ��1֡
        frameout=x(i,:)';            % ����Ҫ�ص����,�����ϳ�����
    else
        M=length(frameout);             % �����Ա����ص���Ӵ���ϳ�����
        frameout=[frameout(1:M-overlap); frameout(M-overlap+1:M).*tempr2+xn(1:overlap).*tempr1; xn(overlap+1:nwin)];
    end
end

