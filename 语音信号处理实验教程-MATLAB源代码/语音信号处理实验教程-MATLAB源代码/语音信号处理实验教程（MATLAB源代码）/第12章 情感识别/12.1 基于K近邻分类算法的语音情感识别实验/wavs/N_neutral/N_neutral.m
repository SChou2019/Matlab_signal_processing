%% ��ȡ���ļ����µ������ļ�������
clc 
clear all
close all

wavefilename = '*.wav';
dr = dir(wavefilename);
neutralVec=zeros(140,length(dr));
for i = 1:length( dr )
    disp(dr(i).name);
    neutralVec(:,i)=featvector(dr(i).name);%���ú���featvector��ȡ����
end
disp(length(dr));
for i=1:size(neutralVec,1)
    neutralVec(i,:)=mapzo(neutralVec(i,:));   %���ú���mapzo���й�һ��
end
save N_neutral neutralVec;