function sumpnn=pnn(p_train,t_train,p_test,t_test)
%pnn:����������
%���������
% p_train��ѵ����������
% t_train��ѵ����������ǩ
% p_test��������������
% t_test��������������ǩ
% ���������
% sumpnn�����ֻ������ʶ����
%% ���������ת��Ϊ����
t_train_temp=t_train;
t_train=ind2vec(t_train);

%% ʹ��newpnn��������PNN SPREADѡȡΪ1.5
Spread=1.5;
net=newpnn(p_train,t_train,Spread)

%% ����Ԥ��δ֪����Ч��
Y2=sim(net,p_test);
Y2c=vec2ind(Y2)
sumpnn=[0 0 0 0 0]; %ÿ�������ȷʶ�����
%ͳ��ʶ����ȷ������ 
for i=1:20
    if Y2c(i)==1
        sumpnn(1)=sumpnn(1)+1;
    end
end
for i=21:40
    if Y2c(i)==2
        sumpnn(2)=sumpnn(2)+1;
    end
end
for i=41:60
    if Y2c(i)==3
        sumpnn(3)=sumpnn(3)+1;
    end
end
for i=61:80
    if Y2c(i)==4
        sumpnn(4)=sumpnn(4)+1;
    end
end
for i=81:100
    if Y2c(i)==5
        sumpnn(5)=sumpnn(5)+1;
    end
end
sumpnn=sumpnn./20;