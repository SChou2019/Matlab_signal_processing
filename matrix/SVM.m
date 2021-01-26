%SVM

clear;
clc;
load("smtrain15.mat");
load("smtest15.mat");
train_data = smtrain15(:,:)
test_data = smtest15(:,:)
train_label = smtrain15(:,1)
test_label = smtest15(:,1)

model = svmtarin(train_label,train_data,"-c100-g 0.01")

[predict_label,accuracy] = svmpredict(train_label,train_data，model)

[predict_label,accuracy] = svmpredict(test_label,test_data，model)
