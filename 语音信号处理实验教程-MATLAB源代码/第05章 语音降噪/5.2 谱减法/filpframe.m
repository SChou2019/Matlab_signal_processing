function frameout=filpframe(x,win,inc)

[nf,len]=size(x);
nx=(nf-1) *inc+len;                 %ԭ�źų���
frameout=zeros(nx,1);
nwin=length(win);                   % ȡ����
if (nwin ~= 1)                           % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
    winx=repmat(win',nf,1);
    x=x./winx;                          % ��ȥ�Ӵ���Ӱ��
    x(find(isinf(x)))=0;                %ȥ����0�õ���Inf
end

sig=zeros((nf-1)*inc+len,1);
for i=1:nf
    start=(i-1)*inc+1;    
    xn=x(i,:)';
    sig(start:start+len-1)=sig(start:start+len-1)+xn;
end
frameout=sig;
