function f=div_frame(input_voice,len,inc)
%������f=fra_div��len,inc,input_voice��
%����˵����len��ÿ֡�ĳ���
%          inc:֡�ƶ�       inc<=len      ��len==inc�����޽���
%          input_voice  :�����ļ����ݣ�Ϊn*1ά����
%�������ܣ���x���ݷֳ�ÿ֡����len�����δ��ƶ�Ϊinc�����ݣ�����fΪ�������ݣ�
%          ÿһ��Ϊһ֡.

input_voice=input_voice(:);%ת��Ϊ������
fh=fix(((size(input_voice,1)-len)/inc)+1);
f=zeros(fh+1,len);
i=1;n=1;
while i<=fh
    j=1;
    while j<=len
        f(i,j)=input_voice(n);
        j=j+1;
        n=n+1;
    end
    n=n-len+inc;
    i=i+1;
end

