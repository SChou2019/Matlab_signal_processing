%% ��ȡ���ļ����µ������ļ�������
clc 
clear all
close all

wavefilename = '*.wav';
dr = dir(wavefilename);
angerVec=zeros(140,length(dr));
for i = 1:length( dr )
    disp(dr(i).name);
    angerVec(:,i)=featvector(dr(i).name);%���ú���featvector��ȡ����
end
disp(length(dr));
for i=1:size(angerVec,1)
    angerVec(i,:)=mapzo(angerVec(i,:));   %���ú���mapzo���й�һ��
end
save W_anger angerVec;