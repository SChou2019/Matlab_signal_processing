function NormMatrix = CMN(Matrix)
%��һ������
[r,c]=size(Matrix);
NormMatrix=zeros(r,c);
for i=1:c
    MatMean=mean(Matrix(:,i));  
    NormMatrix(:,i)=Matrix(:,i)-MatMean;
end
