%�������������ص��׷���ȡ������Ϣ
%˵�����������x��Ƕ��������Ϣ���ź�
%NΪ�ֶγ���
%m0��������ϢΪ0ʱ���ӳ٣�m1��������ϢΪ1ʱ���ӳ�
%len��������Ϣ�ĳ���
%���message����ȡ����������Ϣ
function [ message ] = hide_EchoExtract( x_embeded,N,m0,m1,len)
message = zeros(1,len);
for i=1:len
    x = x_embeded(((i-1)*N+1):(i*N));
    xwhat = rceps(x);
    if(xwhat(m0+1)>xwhat(m1+1))   
        message(1,i) =0;
    else message(1,i)=1;
    end
end

