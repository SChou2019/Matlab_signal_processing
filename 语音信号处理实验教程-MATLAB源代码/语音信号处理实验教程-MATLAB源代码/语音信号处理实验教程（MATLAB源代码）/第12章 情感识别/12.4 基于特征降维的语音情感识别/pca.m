function [trainpca, testpca] = pca(trainsample, test, ReducedDim)
%trainsampleѵ����������
%test������������
%ReducedDim����ά������ά��
[pc,score,latent,tsquare] = princomp(trainsample); %PCA��ά
tranMatrix = pc(:,1:ReducedDim); %ת�þ���
trainpca=trainsample*tranMatrix;%��ά��ѵ����������
testpca=test*tranMatrix;%��ά�������������