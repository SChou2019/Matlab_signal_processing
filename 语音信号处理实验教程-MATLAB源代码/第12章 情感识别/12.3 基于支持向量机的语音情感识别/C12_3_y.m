%ʵ��Ҫ�󣺻���SVM���������ʶ��
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
train(1:30,:)=sampleang(1:30,:); %ÿ����ʮ��������Ϊѵ������
test(1:20,:)=sampleang(31:50,:);%ÿ���ʮ��������Ϊ��������
train(31:60,:)=samplehap(1:30,:);
test(21:40,:)=samplehap(31:50,:);%
train(61:90,:)=sampleneu(1:30,:);
test(41:60,:)=sampleneu(31:50,:);%
train(91:120,:)=samplesad(1:30,:);
test(61:80,:)=samplesad(31:50,:);%
train(121:150,:)=samplefear(1:30,:);
test(81:100,:)=samplefear(31:50,:);%
rate=svmclassfiction(train,test);%����SVM���ຯ��
figure(1)
bar(rate,0.5);
set(gca,'XTickLabel',{'����','����','����','����','����'});
ylabel('ʶ����');
xlabel('���ֻ������');



