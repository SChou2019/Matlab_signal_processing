% �����ʶ������
load rec_data.mat;
fprintf('��ʼʶ��\n');
j = 9;
rec_sph=rdata{j}{1}; % ���ѡ��һ����ʶ��������9��
fprintf('����������ʵֵΪ%d\n',j);
rec_fea = mfcc(rec_sph);  % ������ȡ

% �����ǰ�������ڸ�����hmm��p(X|M)
for i=1:10
  pxsm(i) = viterbi(hmm{i}, rec_fea); 
end
[d,n] = max(pxsm); % �о����������ֵ��Ӧ�������Ϊʶ����
fprintf('������ʶ����Ϊ%d\n',n)
