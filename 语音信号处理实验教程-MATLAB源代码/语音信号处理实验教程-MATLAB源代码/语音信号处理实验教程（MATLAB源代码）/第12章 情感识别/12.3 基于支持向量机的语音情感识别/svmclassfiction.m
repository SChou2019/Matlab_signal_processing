function rate=svmclassfiction(samples,test)   %����������з�����
train1=samples(1:60,:);%������������-���˷���ģ��ѵ������
train2=[samples(1:30,:);samples(61:90,:)];%������������-���Է���ģ��ѵ������
train3=[samples(1:30,:);samples(91:120,:)];%������������-���˷���ģ��ѵ������
train4=[samples(1:30,:);samples(121:150,:)];%������������-���·���ģ��ѵ������
train5=[samples(31:60,:);samples(61:90,:)];%�����������-���Է���ģ��ѵ������
train6=[samples(31:60,:);samples(91:120,:)];%�����������-���˷���ģ��ѵ������
train7=[samples(31:60,:);samples(121:150,:)];%�����������-���·���ģ��ѵ������
train8=[samples(61:90,:);samples(91:120,:)];%������������-���˷���ģ��ѵ������
train9=[samples(61:90,:);samples(121:150,:)];%������������-���·���ģ��ѵ������
train10=[samples(91:120,:);samples(121:150,:)];%�������챯��-���·���ģ��ѵ������
for i=1:30                %�����������
    trainlabel(i)=1;
end
for i=30:60               %�����������
    trainlabel(i)=-1;
end
trainlabel=trainlabel';
svmStruct(1)= svmtrain(train1,trainlabel);    %��������SVM����ģ��
svmStruct(2)= svmtrain(train2,trainlabel);    
svmStruct(3)= svmtrain(train3,trainlabel);   
svmStruct(4)= svmtrain(train4,trainlabel);    
svmStruct(5)= svmtrain(train5,trainlabel);    
svmStruct(6)= svmtrain(train6,trainlabel);    
svmStruct(7)= svmtrain(train7,trainlabel);   
svmStruct(8)= svmtrain(train8,trainlabel);    
svmStruct(9)= svmtrain(train9,trainlabel);    
svmStruct(10)= svmtrain(train10,trainlabel);  
sumang=0; %���������ȷʶ�����
sumhap=0; %���������ȷʶ�����
sumneu=0; %���������ȷʶ�����
sumsad=0; %���������ȷʶ�����
sumfea=0; %���������ȷʶ�����
for i=1:100
    for k=1:5
        votes(k)=0;   %���SVM�����������������涨Ϊĳһ������
    end
    for j=1:10
       C(j) = svmclassify(svmStruct(j),test(i,:));
    end
    if(C(1)==1)    %��һ���о������
         votes(1)=votes(1)+1;  %������л��Ʊ��
    elseif(C(1)==-1)
         votes(2)=votes(2)+1;  %������л��Ʊ��
    end
    if(C(2)==1)    %�ڶ����о������
         votes(1)=votes(1)+1;  %������л��Ʊ��
    elseif(C(2)==-1)
         votes(3)=votes(3)+1;  %������л��Ʊ��
    end
    if(C(3)==1)    %�������о������
         votes(1)=votes(1)+1;  %������л��Ʊ��
    elseif(C(3)==-1)
         votes(4)=votes(4)+1;  %������л��Ʊ��
    end
     if(C(4)==1)    %���ĸ��о������
         votes(1)=votes(1)+1;  %������л��Ʊ��
    elseif(C(4)==-1)
         votes(5)=votes(5)+1;  %������л��Ʊ��
     end
     if(C(5)==1)    %������о������
         votes(2)=votes(2)+1;  %������л��Ʊ��
    elseif(C(5)==-1)
         votes(3)=votes(3)+1;  %������л��Ʊ��
     end
     if(C(6)==1)    %�������о������
         votes(2)=votes(2)+1;  %������л��Ʊ��
    elseif(C(6)==-1)
         votes(4)=votes(4)+1;  %������л��Ʊ��
     end
     if(C(7)==1)    %���߸��о������
         votes(2)=votes(2)+1;  %������л��Ʊ��
    elseif(C(7)==-1)
         votes(5)=votes(5)+1;  %������л��Ʊ��
     end
     if(C(8)==1)    %�ڰ˸��о������
         votes(3)=votes(3)+1;  %������л��Ʊ��
     elseif(C(8)==-1)
         votes(4)=votes(4)+1;  %������л��Ʊ��
     end
     if(C(9)==1)    %�ھŸ��о������
         votes(3)=votes(3)+1;  %������л��Ʊ��
     elseif(C(9)==-1)
        votes(5)=votes(5)+1;  %������л��Ʊ��
     end
     if(C(10)==1)    %��ʮ���о������
         votes(4)=votes(4)+1;  %������л��Ʊ��
     elseif(C(10)==-1)
         votes(5)=votes(5)+1;  %������л��Ʊ��
    end
   if(i>=1&&i<=20&&votes(1)>votes(2)&&votes(1)>votes(3)&&votes(1)>votes(4)&&votes(1)>votes(5))
       sumang=sumang+1;  %������������ȷʶ�����
   end
   if(i>=21&&i<=40&&votes(2)>votes(1)&&votes(2)>votes(3)&&votes(2)>votes(4)&&votes(2)>votes(5))
       sumhap=sumhap+1;  %������������ȷʶ�����
   end
    if(i>=41&&i<=60&&votes(3)>votes(2)&&votes(3)>votes(1)&&votes(3)>votes(4)&&votes(3)>votes(5))
       sumneu=sumneu+1;  %������������ȷʶ�����
    end
    if(i>=61&&i<=80&&votes(4)>votes(1)&&votes(4)>votes(2)&&votes(4)>votes(3)&&votes(4)>votes(5))
       sumsad=sumsad+1;  %������������ȷʶ�����
    end
    if(i>=81&&i<=100&&votes(5)>votes(1)&&votes(5)>votes(3)&&votes(5)>votes(4)&&votes(5)>votes(2))
       sumfea=sumfea+1;  %������������ȷʶ�����
    end
end
rate=[sumang/20,sumhap/20,sumneu/20,sumsad/20,sumfea/20];
rate
bar(rate,0.5);
set(gca,'XTickLabel',{'����','����','����','����','����'});
ylabel('ʶ����');
xlabel('���ֻ������');