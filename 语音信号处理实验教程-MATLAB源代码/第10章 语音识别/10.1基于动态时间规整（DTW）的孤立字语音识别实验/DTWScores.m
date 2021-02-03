function AllScores = DTWScores(rMatrix,N)
%��̬ʱ�������DTW��Ѱ����Сʧ��
%���������rMatrixΪ��ǰ����������MFCC��������,NΪÿ��ģ�������ʻ���
%���������AllScoresΪ

%��ʼ��DTW�б����
Scores1 = zeros(1,N);                
Scores2 = zeros(1,N);
Scores3 = zeros(1,N);


%����ģ������
s1 = load('Vectors1.mat');
fMatrixall1 = struct2cell(s1);
s2 = load('Vectors2.mat');
fMatrixall2 = struct2cell(s2);
s3 = load('Vectors3.mat');
fMatrixall3 = struct2cell(s3);


%����DTW
for i = 1:N
    fMatrix1 = fMatrixall1{i,1};
    fMatrix1 = CMN(fMatrix1);
    Scores1(i) = myDTW(fMatrix1,rMatrix);
end

for j = 1:N
    fMatrix2 = fMatrixall2{j,1};
    fMatrix2 = CMN(fMatrix2);
    Scores2(j) = myDTW(fMatrix2,rMatrix);
end

for k= 1:N
    fMatrix3 = fMatrixall3{k,1};
    fMatrix3 = CMN(fMatrix3);
    Scores3(k) = myDTW(fMatrix3,rMatrix);
end


AllScores = [Scores1;Scores2;Scores3];

