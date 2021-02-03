function [cost] = myDTW(F,R)
%���������FΪģ��MFCC��������RΪ��ǰ����MFCC��������
%���������costΪ���ƥ�����

[r1,c1]=size(F);         %ģ���ά��
[r2,c2]=size(R);         %��ǰ����ά��
distance = zeros(r1,r2);
for n=1:r1
    for m=1:r2
        FR=F(n,:)-R(m,:);
        FR=FR.^2;
        distance(n,m)=sqrt(sum(FR))/c1;     %����ŷ�Ͼ���
    end
end

D = zeros(r1+1,r2+1);   
D(1,:) = inf;        
D(:,1) = inf;          
D(1,1) = 0;
D(2:(r1+1), 2:(r2+1)) = distance;


%Ѱ���������̵����ƥ�����
for i = 1:r1; 
 for j = 1:r2;
   [dmin] = min([D(i, j), D(i, j+1), D(i+1, j)]);
   D(i+1,j+1) = D(i+1,j+1)+dmin;
 end
end

cost = D(r1+1,r2+1);    %�������վ���

