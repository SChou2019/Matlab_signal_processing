clear all;
% ����ѵ�����ݼ�tra_data.mat
load tra_data.mat;

N = 4;   % hmm��״̬��
M = [3,3,3,3]; % ÿ��״̬��Ӧ�Ļ��ģ�ͳɷ���

for i = 1:length(tdata)  % ���ֵ�ѭ��
    fprintf('\n��������%d��mfcc��������\n',i);
    for k = 1:length(tdata{i})  % ��������ѭ��
      obs(k).sph = tdata{i}{k};  % ����i�ĵ�k������
      obs(k).fea = mfcc(obs(k).sph);  % ��������ȡmfcc��������
    end
    
    fprintf('\nѵ������%d��hmm\n',i);
    hmm_temp=inithmm(obs,N,M); %��ʼ��hmmģ��
    hmm{i}=baum_welch(hmm_temp,obs); %��������hmm�ĸ�����
end
fprintf('\nѵ����ɣ�\n');


