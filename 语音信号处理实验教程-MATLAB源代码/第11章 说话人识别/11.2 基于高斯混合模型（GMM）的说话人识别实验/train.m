clear all;
Spk_num=6; %˵���˸���
Tra_num=5;  %ÿ��˵��������ѵ����������Ŀ

ncentres=16; %��ϳɷ���Ŀ
fs=16000; %����Ƶ��

% -- ѵ�� ---
load tra_data.mat; 
for spk_cyc=1:Spk_num    % ����˵����
  fprintf('ѵ����%i��˵����\n',spk_cyc);
  tag1=1;tag2=1; %���ڻ��ܴ洢mfcc
  for sph_cyc=1:Tra_num  % ��������
     speech = tdata{spk_cyc}{sph_cyc}; 
      %---Ԥ����,������ȡ--
     pre_sph=filter([1 -0.97],1,speech); % Ԥ����
     win_type='M'; %������
     cof_num=20; %����ϵ������
     frm_len=fs*0.02; %֡����20ms
     fil_num=20; %�˲��������
     frm_off=fs*0.01; %֡�ƣ�10ms
     c=melcepst(pre_sph,fs,win_type,cof_num,fil_num,frm_len,frm_off); % mfcc������ȡ
     cc=c(:,1:end-1)';
     tag2=tag1+size(cc,2);
     cof(:,tag1:tag2-1)=cc;
     tag1=tag2;
  end
   
  %--- ѵ��GMMģ��---
  kiter=5; %Kmeans������������
  emiter=30; %EM�㷨������������
  mix=gmm_init(ncentres,cof',kiter,'full'); % GMM�ĳ�ʼ��
  [mix,post,errlog]=gmm_em(mix,cof',emiter); % GMM�Ĳ�������
  speaker{spk_cyc}.pai=mix.priors;
  speaker{spk_cyc}.mu=mix.centres;
  speaker{spk_cyc}.sigma=mix.covars;

  clear cof mix;
end
fprintf('ѵ����ɣ�\n',spk_cyc);
save speaker.mat speaker;











