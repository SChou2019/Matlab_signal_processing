%% ��ȡ���ļ����µ������ļ�������
clc 
clear all
close all

wavefilename = '*.wav';
dr = dir(wavefilename);
hapVec=zeros(140,length(dr));
for i = 1:length( dr )
    disp(i)
    disp(dr(i).name);
    hapVec(:,i)=featvector(dr(i).name);%���ú���featvector��ȡ����
end
disp(length(dr));
for i=1:size(hapVec,1)
    hapVec(i,:)=mapzo(hapVec(i,:));   %���ú���mapzo���й�һ��
end
save F_happiness hapVec;