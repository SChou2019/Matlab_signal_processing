%�����ص���ӷ����źŻ�ԭ����
function frameout=Filpframe_OverlapA(x,win,inc)

[nf,len]=size(x);
nx=(nf-1) *inc+len;                 %ԭ�źų���
frameout=zeros(nx,1);
nwin=length(win);                   % ȡ����
if (nwin ~= 1)                           % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
    winx=repmat(win',nf,1);
    x=x./winx;                          % ��ȥ�Ӵ���Ӱ��
    x(find(isinf(x)))=0;                %ȥ����0�õ���Inf
end

for i=1:nf
    start=(i-1)*inc+1;    
    xn=x(i,:)';
    frameout(start:start+len-1)=frameout(start:start+len-1)+xn;
end
