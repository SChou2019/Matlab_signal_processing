%% ��ȡ���ļ����µ������ļ�������
clc 
clear all
close all

wavefilename = '*.wav';
dr = dir(wavefilename);
sadnessVec=zeros(140,length(dr));
for i = 1:length( dr )
    disp(dr(i).name);
    sadnessVec(:,i)=featvector(dr(i).name);%���ú���featvector��ȡ����
end
disp(length(dr));
for i=1:size(sadnessVec,1)
    sadnessVec(i,:)=mapzo(sadnessVec(i,:));   %���ú���mapzo���й�һ��
end
save T_sadness sadnessVec;