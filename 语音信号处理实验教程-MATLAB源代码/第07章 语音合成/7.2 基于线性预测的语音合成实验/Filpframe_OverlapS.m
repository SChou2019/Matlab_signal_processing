%�����ص��洢�����źŻ�ԭ����
function frameout=Filpframe_OverlapS(x,win,inc)

[nf,len]=size(x);
nx=(nf-1) *inc+len;                 %ԭ�źų���
frameout=zeros(nx,1);
nwin=length(win);                   % ȡ����
if (nwin ~= 1)                           % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
    winx=repmat(win',nf,1);
    x=x./winx;                          % ��ȥ�Ӵ���Ӱ��
    x(find(isinf(x)))=0;                %ȥ����0�õ���Inf
end
frameout(1:len)=x(1,:);
for i=1:nf
    frameout(len+(i-1)*inc+1:len+i*inc)=x(i,len-inc+1:len);
end
