%��������������
%˵�����������x��ԭʼ��Ƶ�ź�
%message��������Ϣ
%NΪ�ֶγ���
%m0��������ϢΪ0ʱ���ӳ٣�m1��������ϢΪ1ʱ���ӳ�
%a��˥����
%���x_embeded�Ǻ���������Ϣ���ź�
function [ x_embeded ] = hide_EchoEmbed( x,message,N,m0,m1,a)
x_embeded = x;
nf = min(length(x)/N,length(message));       %����
for i=1:nf
    if(message(i))
        for j=1:N
            if((i-1)*N+j>m1)     x_embeded((i-1)*N+j) = x((i-1)*N+j)+a*x((i-1)*N+j-m1);
            end
        end
    else
        for j=1:N
            if((i-1)*N+j>m0)     x_embeded((i-1)*N+j) = x((i-1)*N+j)+a*x((i-1)*N+j-m0);
            end
        end
    end
end

