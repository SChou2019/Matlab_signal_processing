%ʵ��Ҫ�󣺻���������ά���������ʶ��
clc;
clear;
load A_fear fearVec;
load F_happiness hapVec;
load N_neutral neutralVec;
load T_sadness sadnessVec;
load W_anger angerVec;
sampleang=angerVec';
samplehap=hapVec';
sampleneu=neutralVec';
samplesad=sadnessVec';
samplefear=fearVec';
trainang=sampleang(1:30,:); %ÿ����ʮ��������Ϊѵ������
test(1:20,:)=sampleang(31:50,:);%ÿ���ʮ��������Ϊ��������
trainhap=samplehap(1:30,:);
test(21:40,:)=samplehap(31:50,:);%
trainneu=sampleneu(1:30,:);
test(41:60,:)=sampleneu(31:50,:);%
trainsad=samplesad(1:30,:);
test(61:80,:)=samplesad(31:50,:);%
trainfear=samplefear(1:30,:);
test(81:100,:)=samplefear(31:50,:);%
%��ȡ150������Ϊѵ��������100������ΪԤ������
trainsample=[trainang;trainhap;trainneu;trainsad;trainfear];   %ѵ������
 for i=1:30
   output(i)=1;
 end
 for i=31:60
   output(i)=2;
 end
 for i=61:90
   output(i)=3;
 end
 for i=91:120
   output(i)=4;
 end
 for i=121:150
   output(i)=5;
 end
 trainlabel=output';   %ѵ���������
for i=1:20
   output1(i)=1;
 end
 for i=21:40
   output1(i)=2;
 end
 for i=41:60
   output1(i)=3;
 end
 for i=61:80
   output1(i)=4;
 end
 for i=81:100
   output1(i)=5;
 end
 testlabel=output1';%�����������
[trainpca, testpca] = pca(trainsample, test,5);
rate=knn(trainpca, testpca,7);
figure(1)
bar(rate./20,0.5);
set(gca,'XTickLabel',{'����','����','����','����','����'});
ylabel('ʶ����');
xlabel('���ֻ������');

N=[30,30,30,30,30];
trainsample=[trainang;trainhap;trainneu;trainsad;trainfear];      %ѵ������
testsample(1:20,:)=sampleang(31:50,:);%%ÿ���ʮ��������Ϊ��������
testsample(21:40,:)=samplehap(31:50,:);%
testsample(41:60,:)=sampleneu(31:50,:);%
testsample(61:80,:)=samplesad(31:50,:);%
testsample(81:100,:)=samplefear(31:50,:);%
data=trainsample;
[trainlda, testlda] = lda(data, testsample, N, 4);

rate=knn(trainlda, testlda,7);
figure(2)
bar(rate./20,0.5);
set(gca,'XTickLabel',{'����','����','����','����','����'});
ylabel('ʶ����');
xlabel('���ֻ������');
