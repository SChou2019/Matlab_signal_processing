%% ��ȡ���ļ����µ������ļ�������
clc 
clear all
close all

wavefilename = '*.wav';
dr = dir(wavefilename);
fearVec=zeros(140,length(dr));
for i = 1:length( dr )
    disp(dr(i).name);
    fearVec(:,i)=featvector(dr(i).name);%���ú���featvector��ȡ����
end
disp(length(dr));
for i=1:size(fearVec,1)
    fearVec(i,:)=mapzo(fearVec(i,:));   %���ú���mapzo���й�һ��
end
save A_fear fearVec;