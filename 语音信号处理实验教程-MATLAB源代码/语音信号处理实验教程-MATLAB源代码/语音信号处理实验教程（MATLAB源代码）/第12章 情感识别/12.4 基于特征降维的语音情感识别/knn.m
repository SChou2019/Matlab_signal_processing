function rate= knn(trainsample, test,k)
%ѵ������trainsample
%test��������
sum1=0; %���������ȷʶ�����
sum2=0; %���������ȷʶ�����
sum3=0; %���������ȷʶ�����
sum4=0; %���������ȷʶ�����
sum5=0; %���������ȷʶ�����
summ=0; %��ȷʶ�����
for x=1:100 %50����������
   for y=1:150 %100��ѵ������
       c=(test(x,:)-trainsample(y,:)).^2;
      Eudistance(y)=sqrt(sum(c(:))); 
      %ŷ�Ͼ���
   end;
   [increase,index]=sort(Eudistance);%��������index�����洢����ǰ ��Eudistance�е��±�
   votes1=0;
   votes2=0;
   votes3=0;
   votes4=0;
   votes5=0;
   
   %%%ͶƱ
  
   for n=1:k%�����k��ѵ�����ݵ���ͶƱȨ
      if index(1,n)<=30%��һ��
       votes1=votes1+1;
      elseif index(1,n)>30&&index(1,n)<=60%�ڶ���
       votes2=votes2+1;    
      elseif index(1,n)>60&&index(1,n)<=90%������
       votes3=votes3+1;   
      elseif index(1,n)>91&&index(1,n)<=120%������
       votes4=votes4+1; 
      elseif index(1,n)>120&&index(1,n)<=150%������
       votes5=votes5+1; 
      end
   end
   %%ͶƱ���
   if votes1>=votes2&&votes1>=votes3&&votes1>=votes4&&votes1>=votes5  %ʶ��Ϊ����
       result=1;
   elseif  votes2>=votes1&&votes2>=votes3&&votes2>=votes4&&votes2>=votes5 %ʶ��Ϊ����
       result=2;
   elseif  votes3>=votes1&&votes3>=votes2&&votes3>=votes4&&votes3>=votes5 %ʶ��Ϊ����
       result=3;
    elseif  votes4>=votes1&&votes4>=votes2&&votes4>=votes3&&votes4>=votes5 %ʶ��Ϊ����
       result=4;
     elseif  votes5>=votes1&&votes5>=votes2&&votes5>=votes4&&votes5>=votes3 %ʶ��Ϊ����
       result=5;
   end
    if (x<=20&&result==1)%��ȷʶ��Ϊ��������
        sum1=sum1+1;
    elseif(x>20&&x<=40&&result==2)%��ȷʶ��Ϊ���˸���
        sum2=sum2+1;
    elseif(x>40&&x<=60&&result==3)%��ȷʶ��Ϊ���Ը���
        sum3=sum3+1;
     elseif(x>60&&x<=80&&result==4)%��ȷʶ��Ϊ���˸���
        sum4=sum4+1;
     elseif(x>80&&x<=100&&result==5)%��ȷʶ��Ϊ���¸���
        sum5=sum5+1;
    end
  summ=sum1+sum2+sum3+sum4+sum5;
end
rate=[sum1,sum2,sum3,sum4,sum5];
