function [trainlda, testlda] = lda(data, testsample, N, reduced_dim)
% �������
% data��m*n��ԭʼѵ�����ݣ�mΪ����������nΪά��
% testsample����������
% N����������������������data�е����ݶ�Ӧ
% reduced_dim���µ�����ά��
% �������
% trainlda������LDA������ѵ����������
% testlda������LDA�����Ĳ�����������
C=length(N);
dim=size(data',1);% ����ÿ��������data�е���ʼ����ֹ����
pos=zeros(C,2);
for i=1:C
    START=1;
    if i>1
        START=START+sum(N(1:i-1));
    end
    END=sum(N(1:i));
    pos(i,:)=[START END];
end% ÿ��������ֵ
UI=[];
for i=1:C
    if pos(i,1)==pos(i,2)
        % pos(i,1)==pos(i,2)ʱ��mean�������ܹ���
        UI=[UI;data(pos(i,1),:)];
    else
        UI=[UI;mean(data(pos(i,1):pos(i,2),:))];
    end
end
% �����ֵ
U=mean(data);% ���ɢ�Ⱦ���
SB=zeros(dim,dim);
for i=1:C
    SB=SB+N(i)*(UI(i,:)-U)'*(UI(i,:)-U);
end% ����ɢ�Ⱦ���
SW=zeros(dim,dim);
for i=1:C
    for j=pos(i,1):pos(i,2)
        SW=SW+(data(j,:)-UI(i,:))'*(data(j,:)-UI(i,:));
    end
end% �ò��ֿ���Ҫ��Ҳ���Բ�Ҫ
SW=SW/sum(N);
SB=SB/sum(N);% ��������ֵ����������
matrix=pinv(SW)*SB;
[V,D]=eig(matrix);
condition=dim-reduced_dim+1:dim;
V=V(:,condition);% �����µ�����������������ӳ�䵽�¿ռ�
trainlda=data*V;
testlda=testsample*V;