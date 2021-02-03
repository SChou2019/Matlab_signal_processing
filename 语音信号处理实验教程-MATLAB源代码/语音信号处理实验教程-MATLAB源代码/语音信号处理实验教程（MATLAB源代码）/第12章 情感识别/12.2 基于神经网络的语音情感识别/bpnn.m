function sum=bpnn(trainsample,testsample,class)
%���������trainsample��ѵ������,testsample�ǲ�������,class��ʾѵ�������������trainsample�����ݶ�Ӧ
%sum�����ֻ�����е�ʶ����
for i=1:140
    feature(:,i)= trainsample(:,i);
end
%����ֵ��һ��
[input,minI,maxI] = premnmx( feature')  ;

%�����������
s = length( class ) ;
output = zeros( s , 5  ) ;
for i = 1 : s 
   output( i , class( i )  ) = 1 ;
end

%����������
net = newff( minmax(input) , [10 5] , { 'logsig' 'purelin' } , 'traingdx' ) ;   %����ǰ��������

%����ѵ������
net.trainparam.show = 50 ;
net.trainparam.epochs = 150 ;
net.trainparam.goal = 0.1 ;
net.trainParam.lr = 0.05 ;

%��ʼѵ��
net = train( net, input , output' ) ;

%��ȡ��������
for i=1:140
    featuretest(:,i)= testsample(:,i);
end
 c=testsample(:,141);
%�������ݹ�һ��
testInput = tramnmx(featuretest' , minI, maxI ) ;

%����
Y = sim( net , testInput ) 
sum=[0 0 0 0 0]; %ÿ�������ȷʶ�����
%ͳ��ʶ����ȷ������ 
for i=1:20
    if Y(1,i)>Y(2,i)&&Y(1,i)>Y(3,i)&&Y(1,i)>Y(4,i)&&Y(1,i)>Y(5,i)
        sum(1)=sum(1)+1;
    end
end
for i=21:40
    if Y(2,i)>Y(1,i)&&Y(2,i)>Y(3,i)&&Y(2,i)>Y(4,i)&&Y(2,i)>Y(5,i)
         sum(2)=sum(2)+1;
    end
end
for i=41:60
    if Y(3,i)>Y(1,i)&&Y(3,i)>Y(2,i)&&Y(3,i)>Y(4,i)&&Y(3,i)>Y(5,i)
       sum(3)=sum(3)+1;
    end
end
for i=61:80
    if Y(4,i)>Y(1,i)&&Y(4,i)>Y(2,i)&&Y(4,i)>Y(3,i)&&Y(4,i)>Y(5,i)
       sum(4)=sum(4)+1;
    end
end
for i=81:100
    if Y(5,i)>Y(1,i)&&Y(5,i)>Y(2,i)&&Y(5,i)>Y(3,i)&&Y(5,i)>Y(4,i)
       sum(5)=sum(5)+1;
    end
end
sum=sum./20;